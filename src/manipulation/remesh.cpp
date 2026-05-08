#include "bromesh/manipulation/remesh.h"
#include "bromesh/manipulation/poly_mesh.h"
#include "bromesh/manipulation/weld.h"

#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <vector>

namespace bromesh {

// Isotropic remeshing (Botsch & Kobbelt, 2004 — "A Remeshing Approach to
// Multiresolution Modeling"), implemented over PolyMesh's half-edge surgery
// API (splitEdge / collapseEdge / flipEdge). The hot loop touches no hash
// maps — every topology change is a fixed-cost local rewire of pointers
// inside PolyMesh's flat arrays.
//
// Each pass:
//   1. Split every edge longer than 4/3 · L. splitEdge bisects both adjacent
//      faces (or the lone face on a boundary edge), so neighbors never end
//      up with a T-junction.
//   2. Collapse edges shorter than 4/5 · L, skipping any collapse that
//      would produce an edge longer than 4/3 · L (anti-thrash) or that
//      touches the boundary.
//   3. Flip interior edges toward target valence (6 interior / 4 boundary).
//   4. Tangential Laplacian relaxation; boundary vertices are pinned.

namespace {

inline float edgeLength(const PolyMesh& pm, int32_t hi) {
    const auto& he = pm.halfEdges()[hi];
    int32_t v0 = he.origin;
    int32_t v1 = pm.halfEdges()[he.next].origin;
    float a[3], b[3];
    pm.getVertex(v0, a);
    pm.getVertex(v1, b);
    float dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

// Compute boundary-vertex bitmap in one pass over live half-edges.
std::vector<uint8_t> computeBoundary(const PolyMesh& pm) {
    std::vector<uint8_t> isB(pm.vertexCount(), 0);
    const auto& hes = pm.halfEdges();
    for (size_t i = 0; i < hes.size(); ++i) {
        if (!pm.isLiveHalfEdge((int32_t)i)) continue;
        if (hes[i].twin != PolyMesh::NONE) continue;
        int32_t v0 = hes[i].origin;
        int32_t v1 = hes[hes[i].next].origin;
        if (v0 >= 0 && v0 < (int32_t)isB.size()) isB[v0] = 1;
        if (v1 >= 0 && v1 < (int32_t)isB.size()) isB[v1] = 1;
    }
    return isB;
}

// Snapshot one half-edge per undirected edge. Convention: for an edge with
// twin t, keep min(hi, t); for boundary edges keep hi.
std::vector<int32_t> snapshotEdges(const PolyMesh& pm) {
    std::vector<int32_t> out;
    out.reserve(pm.halfEdgeCount() / 2 + 1);
    const auto& hes = pm.halfEdges();
    for (size_t i = 0; i < hes.size(); ++i) {
        if (!pm.isLiveHalfEdge((int32_t)i)) continue;
        int32_t tw = hes[i].twin;
        if (tw != PolyMesh::NONE && (int32_t)i > tw) continue;
        out.push_back((int32_t)i);
    }
    return out;
}

void splitLongEdges(PolyMesh& pm, float splitThreshold,
                    std::vector<uint8_t>& isBoundaryV) {
    auto edges = snapshotEdges(pm);
    for (int32_t hi : edges) {
        if (!pm.isLiveHalfEdge(hi)) continue;
        if (edgeLength(pm, hi) <= splitThreshold) continue;
        bool wasBoundary = (pm.halfEdges()[hi].twin == PolyMesh::NONE);
        int32_t M = pm.splitEdge(hi);
        if (M < 0) continue;
        if ((int32_t)isBoundaryV.size() <= M)
            isBoundaryV.resize(pm.vertexCount(), 0);
        isBoundaryV[M] = wasBoundary ? 1 : 0;
    }
}

void collapseShortEdges(PolyMesh& pm, float collapseThreshold,
                        float splitThreshold,
                        const std::vector<uint8_t>& isBoundaryV) {
    auto edges = snapshotEdges(pm);
    if (edges.empty()) return;

    // 1-ring fence: never collapse two edges that share a vertex within the
    // same pass. The original implementation's "touchedVert" concept,
    // preserved verbatim — keeps the local topology consistent against
    // overlapping commits.
    std::vector<uint8_t> touched(pm.vertexCount(), 0);

    // Build adjacency once for the anti-thrash check (post-collapse edge
    // length to every neighbor of either endpoint).
    std::vector<std::vector<int32_t>> adj(pm.vertexCount());
    {
        const auto& hes = pm.halfEdges();
        for (size_t i = 0; i < hes.size(); ++i) {
            if (!pm.isLiveHalfEdge((int32_t)i)) continue;
            int32_t a = hes[i].origin;
            int32_t b = hes[hes[i].next].origin;
            auto& la = adj[a];
            bool seen = false;
            for (int32_t x : la) if (x == b) { seen = true; break; }
            if (!seen) la.push_back(b);
        }
    }

    for (int32_t hi : edges) {
        if (!pm.isLiveHalfEdge(hi)) continue;
        if (edgeLength(pm, hi) >= collapseThreshold) continue;

        const auto& he = pm.halfEdges()[hi];
        int32_t a = he.origin;
        int32_t b = pm.halfEdges()[he.next].origin;
        if (a < 0 || b < 0) continue;
        if (touched[a] || touched[b]) continue;
        if (isBoundaryV[a] || isBoundaryV[b]) continue;

        // Anti-thrash: would the merged neighborhood produce an edge
        // longer than the split threshold? If so, splits would just undo
        // the collapse next pass — bail.
        float ax, ay, az, bx, by, bz;
        {
            float p[3]; pm.getVertex(a, p); ax = p[0]; ay = p[1]; az = p[2];
            pm.getVertex(b, p);             bx = p[0]; by = p[1]; bz = p[2];
        }
        const float mx = 0.5f * (ax + bx);
        const float my = 0.5f * (ay + by);
        const float mz = 0.5f * (az + bz);

        bool ok = true;
        auto checkRing = [&](int32_t v) {
            for (int32_t n : adj[v]) {
                if (n == a || n == b) continue;
                float p[3]; pm.getVertex(n, p);
                float dx = p[0] - mx, dy = p[1] - my, dz = p[2] - mz;
                float len = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (len > splitThreshold) { ok = false; return; }
            }
        };
        checkRing(a);
        if (ok) checkRing(b);
        if (!ok) continue;

        const float mid[3] = { mx, my, mz };
        if (!pm.collapseEdge(hi, mid)) continue;

        // Fence the 1-ring of the merged vertex.
        touched[a] = 1;
        touched[b] = 1;
        for (int32_t n : adj[a]) if (n >= 0) touched[n] = 1;
        for (int32_t n : adj[b]) if (n >= 0) touched[n] = 1;
    }
}

int targetValence(int32_t v, const std::vector<uint8_t>& isBoundaryV) {
    return isBoundaryV[v] ? 4 : 6;
}

void flipEdges(PolyMesh& pm, const std::vector<uint8_t>& isBoundaryV) {
    // Per-vertex valence (live outgoing he count).
    std::vector<int> valence(pm.vertexCount(), 0);
    {
        const auto& hes = pm.halfEdges();
        for (size_t i = 0; i < hes.size(); ++i) {
            if (!pm.isLiveHalfEdge((int32_t)i)) continue;
            int32_t v = hes[i].origin;
            if (v >= 0 && v < (int32_t)valence.size()) ++valence[v];
        }
    }

    auto edges = snapshotEdges(pm);
    for (int32_t hi : edges) {
        if (!pm.isLiveHalfEdge(hi)) continue;
        const auto& he = pm.halfEdges()[hi];
        if (he.twin == PolyMesh::NONE) continue;

        // Triangle precondition.
        int32_t h2 = he.next;
        int32_t h3 = pm.halfEdges()[h2].next;
        if (pm.halfEdges()[h3].next != hi) continue;
        int32_t t1 = he.twin;
        int32_t t2 = pm.halfEdges()[t1].next;
        int32_t t3 = pm.halfEdges()[t2].next;
        if (pm.halfEdges()[t3].next != t1) continue;

        int32_t a  = he.origin;
        int32_t b  = pm.halfEdges()[h2].origin;
        int32_t c  = pm.halfEdges()[h3].origin;
        int32_t d  = pm.halfEdges()[t3].origin;
        if (a == d || b == c) continue;
        if (isBoundaryV[a] && isBoundaryV[b]) continue;

        int devBefore = std::abs(valence[a]      - targetValence(a, isBoundaryV))
                      + std::abs(valence[b]      - targetValence(b, isBoundaryV))
                      + std::abs(valence[c]      - targetValence(c, isBoundaryV))
                      + std::abs(valence[d]      - targetValence(d, isBoundaryV));
        int devAfter  = std::abs(valence[a]  - 1 - targetValence(a, isBoundaryV))
                      + std::abs(valence[b]  - 1 - targetValence(b, isBoundaryV))
                      + std::abs(valence[c]  + 1 - targetValence(c, isBoundaryV))
                      + std::abs(valence[d]  + 1 - targetValence(d, isBoundaryV));
        if (devAfter >= devBefore) continue;

        if (!pm.flipEdge(hi)) continue;
        --valence[a]; --valence[b];
        ++valence[c]; ++valence[d];
    }
}

void relax(PolyMesh& pm, const std::vector<uint8_t>& isBoundaryV) {
    const size_t V = pm.vertexCount();
    std::vector<std::vector<int32_t>> adj(V);
    {
        const auto& hes = pm.halfEdges();
        for (size_t i = 0; i < hes.size(); ++i) {
            if (!pm.isLiveHalfEdge((int32_t)i)) continue;
            int32_t a = hes[i].origin;
            int32_t b = hes[hes[i].next].origin;
            auto& la = adj[a];
            bool seen = false;
            for (int32_t x : la) if (x == b) { seen = true; break; }
            if (!seen) la.push_back(b);
        }
    }

    // Per-vertex normal: average of incident face normals, normalized.
    std::vector<float> normals(V * 3, 0.0f);
    for (size_t fi = 0; fi < pm.faceCount(); ++fi) {
        if (!pm.isLiveFace((int32_t)fi)) continue;
        auto fv = pm.faceVertices((int32_t)fi);
        if (fv.size() != 3) continue;
        float p0[3], p1[3], p2[3];
        pm.getVertex(fv[0], p0);
        pm.getVertex(fv[1], p1);
        pm.getVertex(fv[2], p2);
        float e1[3] = { p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2] };
        float e2[3] = { p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2] };
        float nx = e1[1]*e2[2] - e1[2]*e2[1];
        float ny = e1[2]*e2[0] - e1[0]*e2[2];
        float nz = e1[0]*e2[1] - e1[1]*e2[0];
        for (int k = 0; k < 3; ++k) {
            normals[fv[k]*3+0] += nx;
            normals[fv[k]*3+1] += ny;
            normals[fv[k]*3+2] += nz;
        }
    }
    for (size_t v = 0; v < V; ++v) {
        float* n = &normals[v*3];
        float L = std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        if (L > 1e-8f) { n[0]/=L; n[1]/=L; n[2]/=L; }
    }

    // Double-buffer position writes so neighbor updates don't leak.
    std::vector<float> newPos(V * 3);
    for (size_t v = 0; v < V; ++v) {
        float p[3]; pm.getVertex((int32_t)v, p);
        newPos[v*3+0] = p[0];
        newPos[v*3+1] = p[1];
        newPos[v*3+2] = p[2];
    }

    for (size_t v = 0; v < V; ++v) {
        if (adj[v].empty()) continue;
        if (isBoundaryV[v]) continue;

        float cx = 0, cy = 0, cz = 0;
        for (int32_t n : adj[v]) {
            float p[3]; pm.getVertex(n, p);
            cx += p[0]; cy += p[1]; cz += p[2];
        }
        float invN = 1.0f / (float)adj[v].size();
        cx *= invN; cy *= invN; cz *= invN;

        float p[3]; pm.getVertex((int32_t)v, p);
        float dx = cx - p[0];
        float dy = cy - p[1];
        float dz = cz - p[2];

        const float* n = &normals[v*3];
        float dot = dx*n[0] + dy*n[1] + dz*n[2];
        dx -= dot * n[0];
        dy -= dot * n[1];
        dz -= dot * n[2];

        newPos[v*3+0] = p[0] + dx * 0.5f;
        newPos[v*3+1] = p[1] + dy * 0.5f;
        newPos[v*3+2] = p[2] + dz * 0.5f;
    }

    // Apply.
    // PolyMesh's vertices_ is private; we set positions via translateVertex
    // (origin → new pos).
    for (size_t v = 0; v < V; ++v) {
        float cur[3]; pm.getVertex((int32_t)v, cur);
        float off[3] = {
            newPos[v*3+0] - cur[0],
            newPos[v*3+1] - cur[1],
            newPos[v*3+2] - cur[2]
        };
        pm.translateVertex((int32_t)v, off);
    }
}

MeshData polyMeshToTriangles(const PolyMesh& pm) {
    MeshData out;
    out.positions.reserve(pm.vertexCount() * 3);
    for (size_t v = 0; v < pm.vertexCount(); ++v) {
        float p[3]; pm.getVertex((int32_t)v, p);
        out.positions.push_back(p[0]);
        out.positions.push_back(p[1]);
        out.positions.push_back(p[2]);
    }
    out.indices.reserve(pm.faceCount() * 3);
    for (size_t f = 0; f < pm.faceCount(); ++f) {
        if (!pm.isLiveFace((int32_t)f)) continue;
        auto fv = pm.faceVertices((int32_t)f);
        if (fv.size() != 3) continue;
        out.indices.push_back((uint32_t)fv[0]);
        out.indices.push_back((uint32_t)fv[1]);
        out.indices.push_back((uint32_t)fv[2]);
    }
    return out;
}

void recomputeVertexNormals(MeshData& m) {
    const size_t V = m.vertexCount();
    m.normals.assign(V * 3, 0.0f);
    for (size_t t = 0; t < m.triangleCount(); ++t) {
        uint32_t i0 = m.indices[t*3+0];
        uint32_t i1 = m.indices[t*3+1];
        uint32_t i2 = m.indices[t*3+2];
        float e1x = m.positions[i1*3+0] - m.positions[i0*3+0];
        float e1y = m.positions[i1*3+1] - m.positions[i0*3+1];
        float e1z = m.positions[i1*3+2] - m.positions[i0*3+2];
        float e2x = m.positions[i2*3+0] - m.positions[i0*3+0];
        float e2y = m.positions[i2*3+1] - m.positions[i0*3+1];
        float e2z = m.positions[i2*3+2] - m.positions[i0*3+2];
        float nx = e1y*e2z - e1z*e2y;
        float ny = e1z*e2x - e1x*e2z;
        float nz = e1x*e2y - e1y*e2x;
        for (uint32_t idx : {i0, i1, i2}) {
            m.normals[idx*3+0] += nx;
            m.normals[idx*3+1] += ny;
            m.normals[idx*3+2] += nz;
        }
    }
    for (size_t v = 0; v < V; ++v) {
        float* n = &m.normals[v*3];
        float L = std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        if (L > 1e-8f) { n[0]/=L; n[1]/=L; n[2]/=L; }
    }
}

} // anon namespace

MeshData remeshIsotropic(const MeshData& mesh, float targetEdgeLength,
                         int iterations) {
    if (mesh.empty() || mesh.indices.empty() || iterations <= 0) return mesh;

    // Weld first: PolyMesh's surgery operates on vertex-identity manifolds.
    // Position-duplicate vertices (e.g. UV-seam splits in `sphere()`) get
    // their twin pointers paired by PolyMesh::rematchTwins's position pass,
    // but splitEdge / collapseEdge / flipEdge mutate vertex indices in
    // ways that quietly desynchronize seam pairs. Welding collapses each
    // duplicate set to a single vertex up front, which makes the surgery
    // both faster and unambiguous. Attribute streams (UVs, colors,
    // skinning) are dropped by remesh anyway, so the weld is non-lossy
    // here.
    MeshData welded = weldVertices(mesh);
    PolyMesh pm = PolyMesh::fromMeshData(welded.positions, welded.indices);

    if (targetEdgeLength <= 0.0f) {
        double totalLen = 0.0;
        int edgeCount = 0;
        const auto& hes = pm.halfEdges();
        for (size_t i = 0; i < hes.size(); ++i) {
            if (!pm.isLiveHalfEdge((int32_t)i)) continue;
            totalLen += edgeLength(pm, (int32_t)i);
            ++edgeCount;
        }
        if (edgeCount == 0) return welded;
        targetEdgeLength = static_cast<float>(totalLen / edgeCount);
    }
    const float splitThreshold    = targetEdgeLength * (4.0f / 3.0f);
    const float collapseThreshold = targetEdgeLength * (4.0f / 5.0f);

    auto isBoundaryV = computeBoundary(pm);

    for (int it = 0; it < iterations; ++it) {
        splitLongEdges(pm, splitThreshold, isBoundaryV);
        collapseShortEdges(pm, collapseThreshold, splitThreshold, isBoundaryV);
        pm.compact();
        isBoundaryV = computeBoundary(pm);
        flipEdges(pm, isBoundaryV);
        relax(pm, isBoundaryV);
    }

    pm.compact();
    MeshData out = polyMeshToTriangles(pm);
    recomputeVertexNormals(out);
    return out;
}

} // namespace bromesh
