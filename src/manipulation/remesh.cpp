#include "bromesh/manipulation/remesh.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace bromesh {

// Isotropic remeshing (Botsch & Kobbelt, 2004 — "A Remeshing Approach to
// Multiresolution Modeling").
//
// Each pass:
//   1. Split every edge longer than 4/3 · L. Splits are edge-based: a single
//      midpoint is inserted into *every* face sharing the edge, so neighbors
//      never end up with a T-junction that then cracks during relaxation.
//   2. Collapse edges shorter than 4/5 · L, skipping any collapse that would
//      produce an edge longer than 4/3 · L (that edge would just be split
//      again next iteration — the classic split/collapse thrash).
//   3. Flip interior edges toward target valence (6 interior / 4 boundary).
//   4. Tangential Laplacian relaxation with boundary vertices held fixed.

struct HalfEdgeMesh {
    struct Vertex {
        float pos[3];
        bool boundary = false;
    };
    std::vector<Vertex> verts;
    std::vector<std::array<uint32_t, 3>> tris;

    float edgeLength(uint32_t a, uint32_t b) const {
        float dx = verts[a].pos[0] - verts[b].pos[0];
        float dy = verts[a].pos[1] - verts[b].pos[1];
        float dz = verts[a].pos[2] - verts[b].pos[2];
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }
};

struct EdgeKey {
    uint32_t a, b;
    EdgeKey(uint32_t x, uint32_t y)
        : a(std::min(x, y)), b(std::max(x, y)) {}
    bool operator==(const EdgeKey& o) const { return a == o.a && b == o.b; }
};
struct EdgeKeyHash {
    size_t operator()(const EdgeKey& k) const {
        return std::hash<uint64_t>()(
            (static_cast<uint64_t>(k.a) << 32) | k.b);
    }
};

using EdgeFaceMap =
    std::unordered_map<EdgeKey, std::vector<size_t>, EdgeKeyHash>;

static EdgeFaceMap buildEdgeFaces(const HalfEdgeMesh& hm) {
    EdgeFaceMap m;
    m.reserve(hm.tris.size() * 3);
    for (size_t t = 0; t < hm.tris.size(); ++t) {
        const auto& tri = hm.tris[t];
        for (int e = 0; e < 3; ++e) {
            m[EdgeKey(tri[e], tri[(e+1)%3])].push_back(t);
        }
    }
    return m;
}

static void markBoundaryVerts(HalfEdgeMesh& hm) {
    for (auto& v : hm.verts) v.boundary = false;
    EdgeFaceMap m = buildEdgeFaces(hm);
    for (const auto& kv : m) {
        if (kv.second.size() == 1) {
            hm.verts[kv.first.a].boundary = true;
            hm.verts[kv.first.b].boundary = true;
        }
    }
}

// Split every edge longer than splitThreshold, preserving manifoldness by
// inserting each midpoint into every face adjacent to the edge.
static void splitLongEdges(HalfEdgeMesh& hm, float splitThreshold) {
    EdgeFaceMap edgeFaces = buildEdgeFaces(hm);

    auto rmEdge = [&](uint32_t x, uint32_t y, size_t f) {
        auto it = edgeFaces.find(EdgeKey(x, y));
        if (it == edgeFaces.end()) return;
        auto& v = it->second;
        v.erase(std::remove(v.begin(), v.end(), f), v.end());
    };
    auto addEdge = [&](uint32_t x, uint32_t y, size_t f) {
        edgeFaces[EdgeKey(x, y)].push_back(f);
    };

    // Snapshot long edges at start of pass. Splitting never lengthens any
    // edge (midpoints are interior), so candidates stay eligible.
    std::vector<EdgeKey> candidates;
    candidates.reserve(edgeFaces.size());
    for (const auto& kv : edgeFaces) {
        if (hm.edgeLength(kv.first.a, kv.first.b) > splitThreshold)
            candidates.push_back(kv.first);
    }

    for (const EdgeKey& key : candidates) {
        auto it = edgeFaces.find(key);
        if (it == edgeFaces.end() || it->second.empty()) continue;

        uint32_t a = key.a, b = key.b;
        HalfEdgeMesh::Vertex mid;
        for (int c = 0; c < 3; ++c)
            mid.pos[c] = (hm.verts[a].pos[c] + hm.verts[b].pos[c]) * 0.5f;
        // A boundary edge has exactly one adjacent face; its midpoint is
        // still on the boundary.
        mid.boundary = (it->second.size() == 1);
        uint32_t midIdx = static_cast<uint32_t>(hm.verts.size());
        hm.verts.push_back(mid);

        // Copy face list — we'll mutate edgeFaces while iterating.
        std::vector<size_t> faces = it->second;
        for (size_t f : faces) {
            auto& tri = hm.tris[f];
            int ei = -1;
            for (int i = 0; i < 3; ++i) {
                uint32_t x = tri[i], y = tri[(i+1)%3];
                if ((x == a && y == b) || (x == b && y == a)) { ei = i; break; }
            }
            if (ei < 0) continue;
            uint32_t v0 = tri[ei];
            uint32_t v1 = tri[(ei+1)%3];
            uint32_t v2 = tri[(ei+2)%3];

            // Face f: (v0,v1,v2) -> (v0, M, v2).
            // Edges of f that change: (v0,v1) and (v1,v2) leave; (v0,v2) stays.
            rmEdge(v0, v1, f);
            rmEdge(v1, v2, f);
            tri = {v0, midIdx, v2};
            addEdge(v0, midIdx, f);
            addEdge(midIdx, v2, f);

            // New face f2: (M, v1, v2) with winding matching original.
            size_t f2 = hm.tris.size();
            hm.tris.push_back({midIdx, v1, v2});
            addEdge(midIdx, v1, f2);
            addEdge(v1, v2, f2);
            addEdge(v2, midIdx, f2);
        }
    }
}

static void collapseShortEdges(HalfEdgeMesh& hm,
                               float collapseThreshold,
                               float splitThreshold) {
    EdgeFaceMap edgeFaces = buildEdgeFaces(hm);

    std::vector<std::vector<uint32_t>> adj(hm.verts.size());
    auto addUnique = [](std::vector<uint32_t>& list, uint32_t v) {
        for (uint32_t x : list) if (x == v) return;
        list.push_back(v);
    };
    for (const auto& kv : edgeFaces) {
        addUnique(adj[kv.first.a], kv.first.b);
        addUnique(adj[kv.first.b], kv.first.a);
    }

    std::vector<bool> removedVert(hm.verts.size(), false);
    std::vector<bool> touchedVert(hm.verts.size(), false);
    std::unordered_set<size_t> removedTri;

    std::vector<EdgeKey> shortEdges;
    for (const auto& kv : edgeFaces) {
        if (hm.edgeLength(kv.first.a, kv.first.b) < collapseThreshold)
            shortEdges.push_back(kv.first);
    }

    for (const EdgeKey& key : shortEdges) {
        uint32_t a = key.a, b = key.b;
        if (removedVert[a] || removedVert[b]) continue;
        // Non-overlapping commits per pass keep adj/edgeFaces consistent.
        if (touchedVert[a] || touchedVert[b]) continue;
        // Boundary preservation: don't collapse edges that touch boundary
        // vertices. MeshyAI meshes are closed; this costs nothing there.
        if (hm.verts[a].boundary || hm.verts[b].boundary) continue;

        // Link condition: the only vertices adjacent to both a and b must
        // be the 1 or 2 opposite-tip vertices of the faces sharing edge
        // (a,b). Any other common neighbor means collapsing a↔b would
        // create a non-manifold edge or a duplicated triangle.
        auto itEF = edgeFaces.find(EdgeKey(a, b));
        if (itEF == edgeFaces.end()) continue;
        std::unordered_set<uint32_t> tips;
        for (size_t f : itEF->second) {
            const auto& tri = hm.tris[f];
            for (int i = 0; i < 3; ++i)
                if (tri[i] != a && tri[i] != b) tips.insert(tri[i]);
        }
        std::unordered_set<uint32_t> nbA(adj[a].begin(), adj[a].end());
        bool linkOK = true;
        for (uint32_t n : adj[b]) {
            if (n == a) continue;
            if (nbA.count(n) && !tips.count(n)) { linkOK = false; break; }
        }
        if (!linkOK) continue;

        float mx = (hm.verts[a].pos[0] + hm.verts[b].pos[0]) * 0.5f;
        float my = (hm.verts[a].pos[1] + hm.verts[b].pos[1]) * 0.5f;
        float mz = (hm.verts[a].pos[2] + hm.verts[b].pos[2]) * 0.5f;

        // Anti-thrash: if collapsing would produce an edge longer than the
        // split threshold, bail. Otherwise split would just undo us.
        bool ok = true;
        auto checkNeighbors = [&](uint32_t v) {
            for (uint32_t n : adj[v]) {
                if (n == a || n == b) continue;
                if (removedVert[n]) continue;
                float dx = hm.verts[n].pos[0] - mx;
                float dy = hm.verts[n].pos[1] - my;
                float dz = hm.verts[n].pos[2] - mz;
                float len = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (len > splitThreshold) { ok = false; return; }
            }
        };
        checkNeighbors(a);
        if (ok) checkNeighbors(b);
        if (!ok) continue;

        hm.verts[a].pos[0] = mx;
        hm.verts[a].pos[1] = my;
        hm.verts[a].pos[2] = mz;
        removedVert[b] = true;

        for (uint32_t n : adj[b]) {
            if (n == a || n == b) continue;
            addUnique(adj[a], n);
        }

        // Fence off this 1-ring: no adjacent edge collapses this pass.
        touchedVert[a] = true;
        touchedVert[b] = true;
        for (uint32_t n : adj[a]) touchedVert[n] = true;

        for (size_t t = 0; t < hm.tris.size(); ++t) {
            if (removedTri.count(t)) continue;
            auto& tri = hm.tris[t];
            bool touchesA = (tri[0] == a || tri[1] == a || tri[2] == a);
            bool touchesB = (tri[0] == b || tri[1] == b || tri[2] == b);
            if (touchesA && touchesB) {
                removedTri.insert(t);
            } else if (touchesB) {
                for (int i = 0; i < 3; ++i)
                    if (tri[i] == b) tri[i] = a;
            }
        }
    }

    std::vector<std::array<uint32_t, 3>> newTris;
    newTris.reserve(hm.tris.size());
    for (size_t t = 0; t < hm.tris.size(); ++t) {
        if (removedTri.count(t)) continue;
        const auto& tri = hm.tris[t];
        if (tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2]) continue;
        newTris.push_back(tri);
    }
    hm.tris = std::move(newTris);
}

static void flipEdges(HalfEdgeMesh& hm) {
    std::vector<int> valence(hm.verts.size(), 0);
    for (const auto& tri : hm.tris) {
        valence[tri[0]]++;
        valence[tri[1]]++;
        valence[tri[2]]++;
    }

    EdgeFaceMap edgeFaces = buildEdgeFaces(hm);
    auto tgt = [&](uint32_t v) { return hm.verts[v].boundary ? 4 : 6; };

    auto rmEdge = [&](uint32_t x, uint32_t y, size_t f) {
        auto it = edgeFaces.find(EdgeKey(x, y));
        if (it == edgeFaces.end()) return;
        auto& v = it->second;
        v.erase(std::remove(v.begin(), v.end(), f), v.end());
    };
    auto addEdge = [&](uint32_t x, uint32_t y, size_t f) {
        edgeFaces[EdgeKey(x, y)].push_back(f);
    };
    auto containsEdge = [&](size_t f, uint32_t a, uint32_t b) {
        const auto& tri = hm.tris[f];
        bool hasA = false, hasB = false;
        for (int i = 0; i < 3; ++i) {
            if (tri[i] == a) hasA = true;
            if (tri[i] == b) hasB = true;
        }
        return hasA && hasB;
    };

    // Iterate a snapshot of candidate edges so we can mutate edgeFaces live.
    std::vector<EdgeKey> candidates;
    candidates.reserve(edgeFaces.size());
    for (const auto& kv : edgeFaces) candidates.push_back(kv.first);

    // edgeDir: +1 if face f traverses edge as (a→b), -1 if (b→a), 0 if neither.
    auto edgeDir = [&](size_t f, uint32_t a, uint32_t b) -> int {
        const auto& tri = hm.tris[f];
        for (int i = 0; i < 3; ++i) {
            if (tri[i] == a && tri[(i+1)%3] == b) return +1;
            if (tri[i] == b && tri[(i+1)%3] == a) return -1;
        }
        return 0;
    };

    for (const EdgeKey& key : candidates) {
        auto it = edgeFaces.find(key);
        if (it == edgeFaces.end() || it->second.size() != 2) continue;

        uint32_t a = key.a, b = key.b;
        if (hm.verts[a].boundary && hm.verts[b].boundary) continue;

        size_t f0 = it->second[0], f1 = it->second[1];
        if (!containsEdge(f0, a, b) || !containsEdge(f1, a, b)) continue;

        // Classify sides so winding is preserved.
        // FA = face traversing a→b (its opposite vert = vA).
        // FB = face traversing b→a (its opposite vert = vB).
        int d0 = edgeDir(f0, a, b), d1 = edgeDir(f1, a, b);
        if (d0 == 0 || d1 == 0 || d0 == d1) continue;
        size_t FA, FB;
        if (d0 > 0) { FA = f0; FB = f1; }
        else        { FA = f1; FB = f0; }

        uint32_t vA = UINT32_MAX, vB = UINT32_MAX;
        for (int i = 0; i < 3; ++i) {
            uint32_t v = hm.tris[FA][i];
            if (v != a && v != b) vA = v;
            v = hm.tris[FB][i];
            if (v != a && v != b) vB = v;
        }
        if (vA == UINT32_MAX || vB == UINT32_MAX || vA == vB) continue;

        // Don't flip if the diagonal (vA,vB) already exists: would become
        // a 3-face non-manifold edge.
        auto itNew = edgeFaces.find(EdgeKey(vA, vB));
        if (itNew != edgeFaces.end() && !itNew->second.empty()) continue;

        int devBefore = std::abs(valence[a]      - tgt(a))
                      + std::abs(valence[b]      - tgt(b))
                      + std::abs(valence[vA]     - tgt(vA))
                      + std::abs(valence[vB]     - tgt(vB));
        int devAfter  = std::abs(valence[a]  - 1 - tgt(a))
                      + std::abs(valence[b]  - 1 - tgt(b))
                      + std::abs(valence[vA] + 1 - tgt(vA))
                      + std::abs(valence[vB] + 1 - tgt(vB));
        if (devAfter >= devBefore) continue;

        // Commit. New triangles preserve CCW winding of the original quad:
        //   FA: (a, vB, vA)   — a-side
        //   FB: (b, vA, vB)   — b-side
        // The shared diagonal (vA,vB) appears as (vB→vA) in FA and (vA→vB)
        // in FB, which is the correct opposite-direction pairing.
        hm.tris[FA] = {a, vB, vA};
        hm.tris[FB] = {b, vA, vB};
        valence[a]--; valence[b]--;
        valence[vA]++; valence[vB]++;

        // edgeFaces incremental update.
        // FA leaves (a,b),(b,vA); joins (a,vB),(vA,vB); stays (vA,a).
        rmEdge(a, b, FA); rmEdge(b, vA, FA);
        addEdge(a, vB, FA); addEdge(vA, vB, FA);
        // FB leaves (a,b),(a,vB); joins (b,vA),(vA,vB); stays (vB,b).
        rmEdge(a, b, FB); rmEdge(a, vB, FB);
        addEdge(b, vA, FB); addEdge(vA, vB, FB);
        auto itAb = edgeFaces.find(EdgeKey(a, b));
        if (itAb != edgeFaces.end() && itAb->second.empty())
            edgeFaces.erase(itAb);
    }
}

static void relax(HalfEdgeMesh& hm) {
    std::vector<std::vector<uint32_t>> adj(hm.verts.size());
    auto addUnique = [](std::vector<uint32_t>& list, uint32_t v) {
        for (uint32_t x : list) if (x == v) return;
        list.push_back(v);
    };
    for (const auto& tri : hm.tris) {
        for (int i = 0; i < 3; ++i) {
            uint32_t a = tri[i], b = tri[(i+1)%3];
            addUnique(adj[a], b);
            addUnique(adj[b], a);
        }
    }

    std::vector<float> normals(hm.verts.size() * 3, 0.0f);
    for (const auto& tri : hm.tris) {
        float e1[3], e2[3];
        for (int c = 0; c < 3; ++c) {
            e1[c] = hm.verts[tri[1]].pos[c] - hm.verts[tri[0]].pos[c];
            e2[c] = hm.verts[tri[2]].pos[c] - hm.verts[tri[0]].pos[c];
        }
        float nx = e1[1]*e2[2] - e1[2]*e2[1];
        float ny = e1[2]*e2[0] - e1[0]*e2[2];
        float nz = e1[0]*e2[1] - e1[1]*e2[0];
        for (int i = 0; i < 3; ++i) {
            normals[tri[i]*3+0] += nx;
            normals[tri[i]*3+1] += ny;
            normals[tri[i]*3+2] += nz;
        }
    }
    for (size_t v = 0; v < hm.verts.size(); ++v) {
        float* n = &normals[v*3];
        float len = std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        if (len > 1e-8f) { n[0]/=len; n[1]/=len; n[2]/=len; }
    }

    // Double-buffer so neighbor updates don't leak into the same pass.
    std::vector<float> newPos(hm.verts.size() * 3);
    for (size_t v = 0; v < hm.verts.size(); ++v) {
        newPos[v*3+0] = hm.verts[v].pos[0];
        newPos[v*3+1] = hm.verts[v].pos[1];
        newPos[v*3+2] = hm.verts[v].pos[2];
    }

    for (size_t v = 0; v < hm.verts.size(); ++v) {
        if (adj[v].empty()) continue;
        if (hm.verts[v].boundary) continue;

        float cx = 0, cy = 0, cz = 0;
        for (uint32_t n : adj[v]) {
            cx += hm.verts[n].pos[0];
            cy += hm.verts[n].pos[1];
            cz += hm.verts[n].pos[2];
        }
        float invN = 1.0f / adj[v].size();
        cx *= invN; cy *= invN; cz *= invN;

        float dx = cx - hm.verts[v].pos[0];
        float dy = cy - hm.verts[v].pos[1];
        float dz = cz - hm.verts[v].pos[2];

        const float* n = &normals[v * 3];
        float dot = dx*n[0] + dy*n[1] + dz*n[2];
        dx -= dot * n[0];
        dy -= dot * n[1];
        dz -= dot * n[2];

        newPos[v*3+0] = hm.verts[v].pos[0] + dx * 0.5f;
        newPos[v*3+1] = hm.verts[v].pos[1] + dy * 0.5f;
        newPos[v*3+2] = hm.verts[v].pos[2] + dz * 0.5f;
    }
    for (size_t v = 0; v < hm.verts.size(); ++v) {
        hm.verts[v].pos[0] = newPos[v*3+0];
        hm.verts[v].pos[1] = newPos[v*3+1];
        hm.verts[v].pos[2] = newPos[v*3+2];
    }
}

MeshData remeshIsotropic(const MeshData& mesh, float targetEdgeLength,
                         int iterations) {
    if (mesh.empty() || mesh.indices.empty() || iterations <= 0) return mesh;

    HalfEdgeMesh hm;
    size_t vertCount = mesh.vertexCount();
    hm.verts.resize(vertCount);
    for (size_t v = 0; v < vertCount; ++v) {
        hm.verts[v].pos[0] = mesh.positions[v * 3 + 0];
        hm.verts[v].pos[1] = mesh.positions[v * 3 + 1];
        hm.verts[v].pos[2] = mesh.positions[v * 3 + 2];
    }

    size_t triCount = mesh.triangleCount();
    hm.tris.resize(triCount);
    for (size_t t = 0; t < triCount; ++t) {
        hm.tris[t] = {mesh.indices[t*3+0], mesh.indices[t*3+1], mesh.indices[t*3+2]};
    }

    if (targetEdgeLength <= 0.0f) {
        double totalLen = 0.0;
        int edgeCount = 0;
        for (const auto& tri : hm.tris) {
            for (int e = 0; e < 3; ++e) {
                totalLen += hm.edgeLength(tri[e], tri[(e+1)%3]);
                edgeCount++;
            }
        }
        targetEdgeLength = static_cast<float>(totalLen / edgeCount);
    }

    const float splitThreshold    = targetEdgeLength * (4.0f / 3.0f);
    const float collapseThreshold = targetEdgeLength * (4.0f / 5.0f);

    markBoundaryVerts(hm);

    for (int iter = 0; iter < iterations; ++iter) {
        splitLongEdges(hm, splitThreshold);
        collapseShortEdges(hm, collapseThreshold, splitThreshold);
        flipEdges(hm);
        relax(hm);
    }

    // Compact: remove unused vertices, build output MeshData.
    std::vector<uint32_t> remap(hm.verts.size(), UINT32_MAX);
    MeshData result;

    for (const auto& tri : hm.tris) {
        for (int i = 0; i < 3; ++i) {
            if (remap[tri[i]] == UINT32_MAX) {
                remap[tri[i]] = static_cast<uint32_t>(result.positions.size() / 3);
                result.positions.push_back(hm.verts[tri[i]].pos[0]);
                result.positions.push_back(hm.verts[tri[i]].pos[1]);
                result.positions.push_back(hm.verts[tri[i]].pos[2]);
            }
        }
        result.indices.push_back(remap[tri[0]]);
        result.indices.push_back(remap[tri[1]]);
        result.indices.push_back(remap[tri[2]]);
    }

    size_t outVertCount = result.vertexCount();
    result.normals.assign(outVertCount * 3, 0.0f);
    for (size_t t = 0; t < result.triangleCount(); ++t) {
        uint32_t i0 = result.indices[t*3+0];
        uint32_t i1 = result.indices[t*3+1];
        uint32_t i2 = result.indices[t*3+2];
        float e1x = result.positions[i1*3+0]-result.positions[i0*3+0];
        float e1y = result.positions[i1*3+1]-result.positions[i0*3+1];
        float e1z = result.positions[i1*3+2]-result.positions[i0*3+2];
        float e2x = result.positions[i2*3+0]-result.positions[i0*3+0];
        float e2y = result.positions[i2*3+1]-result.positions[i0*3+1];
        float e2z = result.positions[i2*3+2]-result.positions[i0*3+2];
        float nx = e1y*e2z - e1z*e2y;
        float ny = e1z*e2x - e1x*e2z;
        float nz = e1x*e2y - e1y*e2x;
        for (uint32_t idx : {i0, i1, i2}) {
            result.normals[idx*3+0] += nx;
            result.normals[idx*3+1] += ny;
            result.normals[idx*3+2] += nz;
        }
    }
    for (size_t v = 0; v < outVertCount; ++v) {
        float* n = &result.normals[v*3];
        float len = std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        if (len > 1e-8f) { n[0]/=len; n[1]/=len; n[2]/=len; }
    }

    return result;
}

} // namespace bromesh
