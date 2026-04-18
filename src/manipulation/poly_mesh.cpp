#include "bromesh/manipulation/poly_mesh.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace bromesh {

namespace {

// Quantization for position-based twin matching. Matches the scene-editor's
// existing EditMesh.posKey so round-trips through toMeshData/fromMeshData
// stay idempotent.
constexpr float POS_QUANT = 1.0e5f;

struct PosKey {
    int64_t x, y, z;
    bool operator==(const PosKey& o) const { return x==o.x && y==o.y && z==o.z; }
};

struct PosKeyHash {
    size_t operator()(const PosKey& k) const {
        // simple mix; collisions only cause slower lookup, not correctness
        uint64_t h = (uint64_t)k.x * 73856093ULL ^
                     (uint64_t)k.y * 19349663ULL ^
                     (uint64_t)k.z * 83492791ULL;
        return (size_t)h;
    }
};

static PosKey keyOf(float x, float y, float z) {
    return {
        (int64_t)std::llround((double)x * (double)POS_QUANT),
        (int64_t)std::llround((double)y * (double)POS_QUANT),
        (int64_t)std::llround((double)z * (double)POS_QUANT),
    };
}

static void cross3(const float a[3], const float b[3], float out[3]) {
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
}
static float dot3(const float a[3], const float b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
static float length3(const float a[3]) {
    return std::sqrt(dot3(a, a));
}
static void normalize3(float v[3]) {
    float L = length3(v);
    if (L > 1.0e-20f) { v[0]/=L; v[1]/=L; v[2]/=L; }
    else              { v[0]=v[1]=v[2]=0.0f; }
}

} // namespace

// ===========================================================================
// Inspection
// ===========================================================================

int PolyMesh::faceVertexCount(int faceIdx) const {
    if (faceIdx < 0 || faceIdx >= (int)faces_.size()) return 0;
    int32_t start = faces_[faceIdx].halfEdge;
    if (start == NONE) return 0;
    int count = 0;
    int32_t cur = start;
    do {
        ++count;
        cur = halfEdges_[cur].next;
        if (count > (int)halfEdges_.size()) return 0; // malformed guard
    } while (cur != start && cur != NONE);
    return count;
}

std::vector<int32_t> PolyMesh::faceVertices(int faceIdx) const {
    std::vector<int32_t> out;
    if (faceIdx < 0 || faceIdx >= (int)faces_.size()) return out;
    int32_t start = faces_[faceIdx].halfEdge;
    if (start == NONE) return out;
    int32_t cur = start;
    size_t guard = 0;
    do {
        out.push_back(halfEdges_[cur].origin);
        cur = halfEdges_[cur].next;
        if (++guard > halfEdges_.size()) { out.clear(); return out; }
    } while (cur != start && cur != NONE);
    return out;
}

std::vector<int32_t> PolyMesh::faceHalfEdges(int faceIdx) const {
    std::vector<int32_t> out;
    if (faceIdx < 0 || faceIdx >= (int)faces_.size()) return out;
    int32_t start = faces_[faceIdx].halfEdge;
    if (start == NONE) return out;
    int32_t cur = start;
    size_t guard = 0;
    do {
        out.push_back(cur);
        cur = halfEdges_[cur].next;
        if (++guard > halfEdges_.size()) { out.clear(); return out; }
    } while (cur != start && cur != NONE);
    return out;
}

void PolyMesh::getVertex(int vi, float out[3]) const {
    if (vi < 0 || vi >= (int)vertices_.size()) {
        out[0] = out[1] = out[2] = 0.0f;
        return;
    }
    const auto& v = vertices_[vi];
    out[0] = v.x; out[1] = v.y; out[2] = v.z;
}

void PolyMesh::computeFaceNormal(int faceIdx, float out[3]) const {
    out[0] = out[1] = out[2] = 0.0f;
    auto verts = faceVertices(faceIdx);
    if (verts.size() < 3) return;
    // Newell's method: robust even for near-degenerate or slightly non-planar
    // faces. Sum cross products of successive edges; normalize.
    for (size_t i = 0; i < verts.size(); ++i) {
        const Vertex& a = vertices_[verts[i]];
        const Vertex& b = vertices_[verts[(i + 1) % verts.size()]];
        out[0] += (a.y - b.y) * (a.z + b.z);
        out[1] += (a.z - b.z) * (a.x + b.x);
        out[2] += (a.x - b.x) * (a.y + b.y);
    }
    normalize3(out);
}

std::vector<int32_t> PolyMesh::facesInGroup(int groupId) const {
    std::vector<int32_t> out;
    for (size_t i = 0; i < faces_.size(); ++i) {
        if (faces_[i].group == groupId) out.push_back((int32_t)i);
    }
    return out;
}

// ===========================================================================
// Construction
// ===========================================================================

PolyMesh PolyMesh::fromMeshData(const std::vector<float>& positions,
                                 const std::vector<uint32_t>& indices,
                                 const std::vector<int32_t>& triToGroup) {
    PolyMesh pm;
    const size_t vcount = positions.size() / 3;
    pm.vertices_.resize(vcount);
    for (size_t i = 0; i < vcount; ++i) {
        pm.vertices_[i].x = positions[i*3 + 0];
        pm.vertices_[i].y = positions[i*3 + 1];
        pm.vertices_[i].z = positions[i*3 + 2];
        pm.vertices_[i].halfEdge = NONE;
    }

    const size_t tcount = indices.size() / 3;
    pm.faces_.reserve(tcount);
    pm.halfEdges_.reserve(tcount * 3);
    for (size_t t = 0; t < tcount; ++t) {
        Face f;
        f.group = (t < triToGroup.size()) ? triToGroup[t] : -1;
        int32_t fi = (int32_t)pm.faces_.size();
        int32_t h0 = (int32_t)pm.halfEdges_.size();
        int32_t h1 = h0 + 1;
        int32_t h2 = h0 + 2;
        pm.halfEdges_.push_back(HalfEdge{(int32_t)indices[t*3 + 0], NONE, h1, fi});
        pm.halfEdges_.push_back(HalfEdge{(int32_t)indices[t*3 + 1], NONE, h2, fi});
        pm.halfEdges_.push_back(HalfEdge{(int32_t)indices[t*3 + 2], NONE, h0, fi});
        f.halfEdge = h0;
        pm.faces_.push_back(f);
        // vertex outgoing-HE seed
        for (int k = 0; k < 3; ++k) {
            int32_t hi = h0 + k;
            int32_t oi = pm.halfEdges_[hi].origin;
            if (pm.vertices_[oi].halfEdge == NONE) {
                pm.vertices_[oi].halfEdge = hi;
            }
        }
    }

    pm.rematchTwins();
    return pm;
}

PolyMesh PolyMesh::fromPolygon(const std::vector<float>& outerXYZ,
                                const float /*normal*/[3],
                                int32_t group) {
    PolyMesh pm;
    const size_t n = outerXYZ.size() / 3;
    if (n < 3) return pm;
    pm.vertices_.resize(n);
    std::vector<int32_t> verts(n);
    for (size_t i = 0; i < n; ++i) {
        pm.vertices_[i] = Vertex{outerXYZ[i*3], outerXYZ[i*3+1], outerXYZ[i*3+2], NONE};
        verts[i] = (int32_t)i;
    }
    pm.addFace(verts, group);
    // No twin matching — a lone face has all boundary half-edges by design.
    return pm;
}

PolyMesh PolyMesh::fromPolygons(const std::vector<float>& positions,
                                 const std::vector<uint32_t>& polyVerts,
                                 const std::vector<uint32_t>& polyOffsets,
                                 const std::vector<int32_t>& faceGroups) {
    PolyMesh pm;
    const size_t vcount = positions.size() / 3;
    pm.vertices_.resize(vcount);
    for (size_t i = 0; i < vcount; ++i) {
        pm.vertices_[i] = Vertex{positions[i*3], positions[i*3+1], positions[i*3+2], NONE};
    }
    const size_t fcount = polyOffsets.size() > 0 ? polyOffsets.size() - 1 : 0;
    for (size_t f = 0; f < fcount; ++f) {
        const uint32_t start = polyOffsets[f];
        const uint32_t end   = polyOffsets[f + 1];
        if (end - start < 3) continue;
        std::vector<int32_t> verts;
        verts.reserve(end - start);
        for (uint32_t k = start; k < end; ++k) verts.push_back((int32_t)polyVerts[k]);
        pm.addFace(verts, f < faceGroups.size() ? faceGroups[f] : -1);
    }
    pm.rematchTwins();
    return pm;
}

// ===========================================================================
// Primitive mutators — addVertex / addFace
// ===========================================================================

int32_t PolyMesh::addVertex(float x, float y, float z) {
    vertices_.push_back(Vertex{x, y, z, NONE});
    return (int32_t)vertices_.size() - 1;
}

int32_t PolyMesh::addFace(const std::vector<int32_t>& verts, int32_t group) {
    if (verts.size() < 3) return NONE;
    const int32_t fi   = (int32_t)faces_.size();
    const int32_t hBase = (int32_t)halfEdges_.size();
    const int32_t N    = (int32_t)verts.size();
    for (int32_t k = 0; k < N; ++k) {
        HalfEdge he;
        he.origin = verts[k];
        he.twin   = NONE;
        he.next   = hBase + ((k + 1) % N);
        he.face   = fi;
        halfEdges_.push_back(he);
        if (vertices_[verts[k]].halfEdge == NONE) {
            vertices_[verts[k]].halfEdge = hBase + k;
        }
    }
    Face f;
    f.halfEdge = hBase;
    f.group    = group;
    faces_.push_back(f);
    return fi;
}

// ===========================================================================
// Twin matching
// ===========================================================================

void PolyMesh::rematchTwins() {
    // Pass 1: vertex-identity pairing.
    // Key = (origin, destination) of each directed half-edge.
    std::unordered_map<uint64_t, int32_t> byIdDir;
    byIdDir.reserve(halfEdges_.size());
    auto idKey = [](int32_t a, int32_t b) -> uint64_t {
        return ((uint64_t)(uint32_t)a << 32) | (uint64_t)(uint32_t)b;
    };
    // Don't clear existing valid twins — just walk unpaired and try.
    for (size_t i = 0; i < halfEdges_.size(); ++i) {
        auto& he = halfEdges_[i];
        if (he.twin != NONE) continue;
        uint64_t k = idKey(he.origin, halfEdges_[he.next].origin);
        byIdDir.emplace(k, (int32_t)i);
    }
    for (size_t i = 0; i < halfEdges_.size(); ++i) {
        auto& he = halfEdges_[i];
        if (he.twin != NONE) continue;
        uint64_t revKey = idKey(halfEdges_[he.next].origin, he.origin);
        auto it = byIdDir.find(revKey);
        if (it == byIdDir.end()) continue;
        int32_t tj = it->second;
        if (tj == (int32_t)i) continue;
        if (halfEdges_[tj].twin != NONE) continue;
        halfEdges_[i].twin = tj;
        halfEdges_[tj].twin = (int32_t)i;
    }

    // Pass 2: quantized-position pairing for hard-edge seams where two
    // faces use distinct vertex indices at the same world point.
    struct PosDirKey { PosKey a, b; };
    struct PosDirKeyHash {
        size_t operator()(const PosDirKey& k) const {
            PosKeyHash h;
            return h(k.a) ^ (h(k.b) * 31ULL);
        }
    };
    struct PosDirKeyEq {
        bool operator()(const PosDirKey& l, const PosDirKey& r) const {
            return l.a == r.a && l.b == r.b;
        }
    };
    std::unordered_map<PosDirKey, int32_t, PosDirKeyHash, PosDirKeyEq> byPosDir;
    auto posKeyOfVert = [&](int32_t vi) {
        const Vertex& v = vertices_[vi];
        return keyOf(v.x, v.y, v.z);
    };
    for (size_t i = 0; i < halfEdges_.size(); ++i) {
        if (halfEdges_[i].twin != NONE) continue;
        PosDirKey k{
            posKeyOfVert(halfEdges_[i].origin),
            posKeyOfVert(halfEdges_[halfEdges_[i].next].origin),
        };
        // First unpaired he with this key wins; duplicates (non-manifold by
        // position) can only pair one twin each.
        byPosDir.emplace(k, (int32_t)i);
    }
    for (size_t i = 0; i < halfEdges_.size(); ++i) {
        if (halfEdges_[i].twin != NONE) continue;
        PosDirKey revK{
            posKeyOfVert(halfEdges_[halfEdges_[i].next].origin),
            posKeyOfVert(halfEdges_[i].origin),
        };
        auto it = byPosDir.find(revK);
        if (it == byPosDir.end()) continue;
        int32_t tj = it->second;
        if (tj == (int32_t)i) continue;
        if (halfEdges_[tj].twin != NONE) continue;
        halfEdges_[i].twin = tj;
        halfEdges_[tj].twin = (int32_t)i;
    }
}

// ===========================================================================
// Tessellation — project each face to its plane, ear-clip in 2D
// ===========================================================================

namespace {

// 2D signed area of triangle a-b-c (positive = CCW).
static float tri2DSignedArea(const float a[2], const float b[2], const float c[2]) {
    return 0.5f * ((b[0]-a[0]) * (c[1]-a[1]) - (b[1]-a[1]) * (c[0]-a[0]));
}

static bool pointInTri2D(const float a[2], const float b[2], const float c[2],
                          const float p[2]) {
    // Barycentric test. Reject on boundary so ear-clipper doesn't pick ears
    // with co-located interior vertices.
    float d1 = tri2DSignedArea(a, b, p);
    float d2 = tri2DSignedArea(b, c, p);
    float d3 = tri2DSignedArea(c, a, p);
    bool hasNeg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    bool hasPos = (d1 > 0) || (d2 > 0) || (d3 > 0);
    return !(hasNeg && hasPos);
}

// Ear-clipping triangulation for a simple CCW polygon in 2D. Writes
// triangles as (i, j, k) index triples referencing uv/idx entries.
// No Steiner points — assumes input polygon is simple and non-self-
// intersecting. Returns false if unable to triangulate (gives up on a
// degenerate/self-intersecting input rather than infinite-looping).
static bool earClip2D(const std::vector<float>& uv,      // stride 2
                      std::vector<uint32_t>& outTris) {
    const int N = (int)(uv.size() / 2);
    if (N < 3) return false;

    // Compute signed area to decide winding.
    float area2 = 0;
    for (int i = 0; i < N; ++i) {
        const float* p = &uv[i*2];
        const float* q = &uv[((i + 1) % N) * 2];
        area2 += p[0]*q[1] - q[0]*p[1];
    }
    const bool ccw = area2 > 0;

    // Working vertex ring as an index list. If CW, walk in reverse so the
    // ear test (requires CCW) works.
    std::vector<int> V(N);
    if (ccw) { for (int i = 0; i < N; ++i) V[i] = i; }
    else     { for (int i = 0; i < N; ++i) V[i] = N - 1 - i; }

    auto at = [&](int idx, float out[2]) {
        out[0] = uv[idx*2 + 0];
        out[1] = uv[idx*2 + 1];
    };

    int remaining = N;
    int guard     = 2 * N;        // degenerate input → bail
    while (remaining > 2) {
        if (--guard < 0) return false;
        bool eared = false;
        for (int i = 0; i < remaining; ++i) {
            int ia = V[(i + remaining - 1) % remaining];
            int ib = V[i];
            int ic = V[(i + 1) % remaining];
            float a[2], b[2], c[2];
            at(ia, a); at(ib, b); at(ic, c);
            // Must be convex (signed area > 0 for CCW).
            if (tri2DSignedArea(a, b, c) <= 0) continue;
            // No other polygon vertex inside.
            bool clean = true;
            for (int j = 0; j < remaining && clean; ++j) {
                int idj = V[j];
                if (idj == ia || idj == ib || idj == ic) continue;
                float p[2]; at(idj, p);
                if (pointInTri2D(a, b, c, p)) clean = false;
            }
            if (!clean) continue;
            outTris.push_back((uint32_t)ia);
            outTris.push_back((uint32_t)ib);
            outTris.push_back((uint32_t)ic);
            V.erase(V.begin() + i);
            --remaining;
            eared = true;
            break;
        }
        if (!eared) return false;   // couldn't find an ear (malformed poly)
    }
    return true;
}

} // namespace

PolyMesh::Tessellation PolyMesh::tessellate() const {
    Tessellation out;
    if (faces_.empty()) return out;

    // Emit vertices per-face to keep per-face flat normals clean (no sharing
    // across face boundaries). Caller can weld if desired.
    out.positions.reserve(vertices_.size() * 3);
    out.normals.reserve(vertices_.size() * 3);
    out.indices.reserve(faces_.size() * 3);
    out.triToFace.reserve(faces_.size());
    out.triToGroup.reserve(faces_.size());

    for (size_t f = 0; f < faces_.size(); ++f) {
        auto verts = faceVertices((int)f);
        if (verts.size() < 3) continue;

        // Pick plane basis from face normal.
        float n[3]; computeFaceNormal((int)f, n);
        if (dot3(n, n) < 1e-20f) continue;

        // Build u, v on the plane. Choose u ⟂ n avoiding the most-aligned axis.
        float ax = std::abs(n[0]), ay = std::abs(n[1]), az = std::abs(n[2]);
        float tmp[3] = {0,0,0};
        if (ax <= ay && ax <= az)      { tmp[0] = 1; }
        else if (ay <= ax && ay <= az) { tmp[1] = 1; }
        else                           { tmp[2] = 1; }
        float u[3]; cross3(tmp, n, u); normalize3(u);
        float v[3]; cross3(n, u, v);   normalize3(v);

        // Project to 2D.
        std::vector<float> uv(verts.size() * 2);
        for (size_t i = 0; i < verts.size(); ++i) {
            const Vertex& vx = vertices_[verts[i]];
            float p[3] = {vx.x, vx.y, vx.z};
            uv[i*2 + 0] = dot3(p, u);
            uv[i*2 + 1] = dot3(p, v);
        }

        std::vector<uint32_t> tris;
        bool ok = earClip2D(uv, tris);
        if (!ok || tris.empty()) {
            // Fallback: naive fan from vertex 0. Still better than dropping
            // the face entirely — callers can detect degenerate output via
            // validate() and computeFaceNormal.
            tris.clear();
            for (size_t i = 1; i + 1 < verts.size(); ++i) {
                tris.push_back(0);
                tris.push_back((uint32_t)i);
                tris.push_back((uint32_t)(i + 1));
            }
        }

        // Emit positions for this face (per-face dup so shading seams work)
        uint32_t baseIdx = (uint32_t)(out.positions.size() / 3);
        for (int32_t vi : verts) {
            const Vertex& vx = vertices_[vi];
            out.positions.push_back(vx.x);
            out.positions.push_back(vx.y);
            out.positions.push_back(vx.z);
            out.normals.push_back(n[0]);
            out.normals.push_back(n[1]);
            out.normals.push_back(n[2]);
        }
        for (size_t t = 0; t < tris.size(); t += 3) {
            out.indices.push_back(baseIdx + tris[t + 0]);
            out.indices.push_back(baseIdx + tris[t + 1]);
            out.indices.push_back(baseIdx + tris[t + 2]);
            out.triToFace.push_back((int32_t)f);
            out.triToGroup.push_back(faces_[f].group);
        }
    }
    return out;
}

// ===========================================================================
// Validation
// ===========================================================================

PolyMesh::Validation PolyMesh::validate() const {
    Validation r;
    for (size_t i = 0; i < halfEdges_.size(); ++i) {
        const HalfEdge& he = halfEdges_[i];
        if (he.origin < 0 || he.origin >= (int32_t)vertices_.size()) {
            r.valid = false;
            r.errors.push_back("he[" + std::to_string(i) + "] has invalid origin");
        }
        if (he.next < 0 || he.next >= (int32_t)halfEdges_.size()) {
            r.valid = false;
            r.errors.push_back("he[" + std::to_string(i) + "] has invalid next");
        }
        if (he.face < 0 || he.face >= (int32_t)faces_.size()) {
            r.valid = false;
            r.errors.push_back("he[" + std::to_string(i) + "] has invalid face");
        }
        if (he.twin == NONE) {
            ++r.boundaryHalfEdges;
        } else if (he.twin < 0 || he.twin >= (int32_t)halfEdges_.size()) {
            r.valid = false;
            r.errors.push_back("he[" + std::to_string(i) + "] has out-of-range twin");
        } else if (halfEdges_[he.twin].twin != (int32_t)i) {
            r.valid = false;
            r.errors.push_back("he[" + std::to_string(i) + "] twin link is not symmetric");
        }
    }
    for (size_t fi = 0; fi < faces_.size(); ++fi) {
        const Face& f = faces_[fi];
        if (f.halfEdge == NONE) {
            r.valid = false;
            r.errors.push_back("face[" + std::to_string(fi) + "] has no halfEdge");
            continue;
        }
        // Walk and check closure + face-pointer consistency.
        int32_t cur = f.halfEdge;
        size_t guard = 0;
        do {
            if (halfEdges_[cur].face != (int32_t)fi) {
                r.valid = false;
                r.errors.push_back("face[" + std::to_string(fi) + "] has he whose face pointer disagrees");
                break;
            }
            cur = halfEdges_[cur].next;
            if (++guard > halfEdges_.size()) {
                r.valid = false;
                r.errors.push_back("face[" + std::to_string(fi) + "] loop doesn't close");
                break;
            }
        } while (cur != f.halfEdge);
    }
    r.isClosed = (r.boundaryHalfEdges == 0);
    return r;
}

// ===========================================================================
// Boundary discovery
// ===========================================================================

std::vector<std::vector<int32_t>> PolyMesh::findFaceBoundary(int faceIdx) const {
    std::vector<std::vector<int32_t>> loops;
    if (faceIdx < 0 || faceIdx >= (int)faces_.size()) return loops;
    // Single face: its boundary *is* the ordered ring of its half-edges
    // whose twin's face is anything other than itself.
    std::vector<int32_t> loop;
    int32_t start = faces_[faceIdx].halfEdge;
    int32_t cur = start;
    size_t guard = 0;
    do {
        loop.push_back(cur);
        cur = halfEdges_[cur].next;
        if (++guard > halfEdges_.size()) { loop.clear(); break; }
    } while (cur != start);
    if (!loop.empty()) loops.push_back(std::move(loop));
    return loops;
}

std::vector<std::vector<int32_t>> PolyMesh::findGroupBoundary(int groupId) const {
    std::vector<std::vector<int32_t>> loops;
    // Seed set: he in `groupId` whose twin is in a different group (or NONE).
    std::vector<int32_t> seeds;
    auto isGroupBoundary = [&](int32_t hi) {
        const HalfEdge& he = halfEdges_[hi];
        if (he.twin == NONE) return true;
        return faces_[halfEdges_[he.twin].face].group != groupId;
    };
    for (size_t i = 0; i < halfEdges_.size(); ++i) {
        if (faces_[halfEdges_[i].face].group != groupId) continue;
        if (isGroupBoundary((int32_t)i)) seeds.push_back((int32_t)i);
    }
    // Walk each seed as the start of a loop; use one-ring via twin.next
    // around the destination vertex to find the next boundary he in group.
    std::unordered_set<int32_t> visited;
    for (int32_t start : seeds) {
        if (visited.count(start)) continue;
        std::vector<int32_t> loop;
        int32_t cur = start;
        while (cur != NONE && !visited.count(cur)) {
            visited.insert(cur);
            loop.push_back(cur);
            // Walk one-ring around he.next.origin in the same group.
            int32_t next = halfEdges_[cur].next;
            int32_t chosen = NONE;
            size_t guard = 0;
            while (next != NONE) {
                if (faces_[halfEdges_[next].face].group == groupId &&
                    isGroupBoundary(next)) { chosen = next; break; }
                int32_t tw = halfEdges_[next].twin;
                if (tw == NONE) { chosen = NONE; break; }
                next = halfEdges_[tw].next;
                if (++guard > halfEdges_.size()) { chosen = NONE; break; }
            }
            cur = chosen;
        }
        if (!loop.empty()) loops.push_back(std::move(loop));
    }
    return loops;
}

// ===========================================================================
// Surgery — translate / extrude
// ===========================================================================

void PolyMesh::translateVertex(int32_t vi, const float offset[3]) {
    if (vi < 0 || vi >= (int32_t)vertices_.size()) return;
    vertices_[vi].x += offset[0];
    vertices_[vi].y += offset[1];
    vertices_[vi].z += offset[2];
}

void PolyMesh::translateFace(int32_t faceIdx, const float offset[3]) {
    auto verts = faceVertices(faceIdx);
    for (int32_t vi : verts) {
        vertices_[vi].x += offset[0];
        vertices_[vi].y += offset[1];
        vertices_[vi].z += offset[2];
    }
}

void PolyMesh::translateFaceWithRing(int32_t faceIdx, const float offset[3]) {
    if (faceIdx < 0 || faceIdx >= (int32_t)faces_.size()) return;
    // For each face whose group matches `faceIdx`'s group — pushes the whole
    // face-group, consistent with scene-editor semantics where a face group
    // is a logical face even if it's split into multiple triangle faces.
    int32_t gIdx = faces_[faceIdx].group;

    // Gather every vertex touching the group, plus any ring-walkable
    // seam-duplicate vert sharing position on a boundary.
    std::unordered_set<int32_t> moved;

    auto collectGroupVerts = [&](std::vector<int32_t>& out) {
        std::unordered_set<int32_t> seen;
        for (size_t fi = 0; fi < faces_.size(); ++fi) {
            if (faces_[fi].group != gIdx) continue;
            auto vs = faceVertices((int)fi);
            for (int32_t v : vs) if (seen.insert(v).second) out.push_back(v);
        }
    };

    // Find group boundary to locate which of the gathered verts are on the
    // boundary (they need ring-walked seam sweeps); interior verts translate
    // once.
    auto boundary = findGroupBoundary(gIdx);
    std::unordered_set<int32_t> boundaryVerts;
    for (const auto& loop : boundary) {
        for (int32_t hi : loop) {
            boundaryVerts.insert(halfEdges_[hi].origin);
            boundaryVerts.insert(halfEdges_[halfEdges_[hi].next].origin);
        }
    }

    std::vector<int32_t> groupVerts;
    collectGroupVerts(groupVerts);

    // Walk a one-ring via twin.next around a boundary vertex, collecting
    // vertices at the same position (seam duplicates).
    auto walkRing = [&](int32_t seedHe, std::vector<int32_t>& ring) {
        std::unordered_set<int32_t> seen;
        int32_t cur = seedHe;
        size_t guard = 0;
        while (cur != NONE && !seen.count(cur)) {
            seen.insert(cur);
            int32_t vo = halfEdges_[cur].origin;
            if (moved.insert(vo).second) ring.push_back(vo);
            int32_t tw = halfEdges_[cur].twin;
            if (tw == NONE) break;
            cur = halfEdges_[tw].next;
            if (++guard > halfEdges_.size()) break;
        }
    };

    std::vector<std::vector<int32_t>> snapGroups;
    // Walk boundary half-edges first so seam dups in adjacent faces are
    // picked up (they won't appear as vertices of the pushed group's faces).
    for (const auto& loop : boundary) {
        for (int32_t hi : loop) {
            std::vector<int32_t> ring;
            walkRing(hi, ring);
            if (ring.size() > 1) snapGroups.push_back(std::move(ring));
            else if (ring.size() == 1) {
                // singleton — still record so it moves
            }
        }
    }
    // Interior verts of the group translate without ring sweep.
    for (int32_t v : groupVerts) {
        if (boundaryVerts.count(v)) continue;   // boundary handled above
        moved.insert(v);
    }

    for (int32_t v : moved) {
        vertices_[v].x += offset[0];
        vertices_[v].y += offset[1];
        vertices_[v].z += offset[2];
    }
    // Snap each ring to one canonical post-move position to avoid Float32
    // seam jitter that would break position-based twin matching on a
    // round-trip through fromMeshData.
    for (const auto& ring : snapGroups) {
        const Vertex& a = vertices_[ring[0]];
        for (size_t i = 1; i < ring.size(); ++i) {
            vertices_[ring[i]].x = a.x;
            vertices_[ring[i]].y = a.y;
            vertices_[ring[i]].z = a.z;
        }
    }
}

PolyMesh::ExtrudeResult PolyMesh::extrudeFace(int32_t faceIdx,
                                              const float offset[3],
                                              bool withBackFace,
                                              int32_t bridgeGroup,
                                              int32_t backGroup) {
    ExtrudeResult out;
    if (faceIdx < 0 || faceIdx >= (int32_t)faces_.size()) return out;

    // Snapshot the ordered boundary of the face BEFORE we mutate half-edges.
    std::vector<int32_t> bdHEs = faceHalfEdges(faceIdx);
    if (bdHEs.size() < 3) return out;

    // Capture each boundary edge's original endpoints and adjacent group
    // before we sever twins.
    struct BoundaryRec {
        int32_t he;
        int32_t oldA, oldB;
        int32_t newA = NONE, newB = NONE;
        int32_t adjGroup = -1;
    };
    std::vector<BoundaryRec> recs;
    recs.reserve(bdHEs.size());
    for (int32_t hi : bdHEs) {
        BoundaryRec r;
        r.he = hi;
        r.oldA = halfEdges_[hi].origin;
        r.oldB = halfEdges_[halfEdges_[hi].next].origin;
        r.adjGroup = (halfEdges_[hi].twin == NONE)
                   ? -1
                   : faces_[halfEdges_[halfEdges_[hi].twin].face].group;
        recs.push_back(r);
    }

    // Duplicate every boundary vertex (append new verts at vert + offset).
    // Use a map from old vertex index to new — boundary vertices may repeat
    // around the loop.
    std::unordered_map<int32_t, int32_t> dupMap;
    auto dupOf = [&](int32_t vi) -> int32_t {
        auto it = dupMap.find(vi);
        if (it != dupMap.end()) return it->second;
        const Vertex& v = vertices_[vi];
        int32_t ni = addVertex(v.x + offset[0], v.y + offset[1], v.z + offset[2]);
        dupMap.emplace(vi, ni);
        return ni;
    };
    for (auto& r : recs) {
        r.newA = dupOf(r.oldA);
        r.newB = dupOf(r.oldB);
    }
    out.dupVerts.reserve(dupMap.size());
    for (auto& kv : dupMap) out.dupVerts.push_back(kv.second);

    // Rewire the face's half-edges to use the duplicate vertices. The face
    // itself now floats one offset away.
    for (int32_t hi : bdHEs) {
        halfEdges_[hi].origin = dupMap.at(halfEdges_[hi].origin);
        if (vertices_[halfEdges_[hi].origin].halfEdge == NONE) {
            vertices_[halfEdges_[hi].origin].halfEdge = hi;
        }
    }

    // Sever twins along the boundary — neighbours (if any) lose their
    // connection to this face until the bridges reconnect.
    for (int32_t hi : bdHEs) {
        int32_t tw = halfEdges_[hi].twin;
        if (tw != NONE) {
            halfEdges_[tw].twin = NONE;
            halfEdges_[hi].twin = NONE;
        }
    }

    // Emit one QUAD bridge face per boundary edge:
    //   quad verts (CCW from the bridge's outward side):
    //     oldA → oldB → newB → newA
    int32_t groupCursor = bridgeGroup;
    out.bridgeFaces.reserve(recs.size());
    out.bridgeAdjGroup.reserve(recs.size());
    for (auto& r : recs) {
        std::vector<int32_t> quad{ r.oldA, r.oldB, r.newB, r.newA };
        int32_t g = (bridgeGroup == -1) ? groupCursor-- : bridgeGroup;
        int32_t fi = addFace(quad, g);
        out.bridgeFaces.push_back(fi);
        out.bridgeAdjGroup.push_back(r.adjGroup);
    }

    if (withBackFace) {
        // Back face = original vertex ring in reverse order, closing the slab.
        std::vector<int32_t> back;
        back.reserve(recs.size());
        for (auto it = recs.rbegin(); it != recs.rend(); ++it) back.push_back(it->oldA);
        out.backFace = addFace(back, backGroup);
    }

    rematchTwins();
    return out;
}

// ===========================================================================
// mergeFacesByGroup — collapse coplanar tri groups into N-gon faces
// ===========================================================================

void PolyMesh::mergeFacesByGroup() {
    if (faces_.empty()) return;

    // Step 1: union-find over faces that share an interior edge AND share a
    // group tag. Two faces F0, F1 are in the same component if there is a
    // half-edge H in F0 with twin H' in F1 and faces_[F0].group == faces_[F1].group.
    std::vector<int32_t> parent(faces_.size());
    for (size_t i = 0; i < parent.size(); ++i) parent[i] = (int32_t)i;
    std::function<int32_t(int32_t)> find = [&](int32_t x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };
    auto unite = [&](int32_t a, int32_t b) {
        int32_t ra = find(a), rb = find(b);
        if (ra != rb) parent[rb] = ra;
    };
    for (size_t i = 0; i < halfEdges_.size(); ++i) {
        const HalfEdge& he = halfEdges_[i];
        if (he.twin == NONE) continue;
        int32_t fa = he.face;
        int32_t fb = halfEdges_[he.twin].face;
        if (fa < 0 || fb < 0) continue;
        if (faces_[fa].group != faces_[fb].group) continue;
        if (faces_[fa].group < 0) continue;   // ungrouped: never merge
        unite(fa, fb);
    }

    // Step 2: collect components.
    std::unordered_map<int32_t, std::vector<int32_t>> compMembers;
    for (size_t i = 0; i < faces_.size(); ++i) {
        compMembers[find((int32_t)i)].push_back((int32_t)i);
    }
    // Skip components with only one face — nothing to merge.
    bool anyMerge = false;
    for (auto& kv : compMembers) if (kv.second.size() > 1) { anyMerge = true; break; }
    if (!anyMerge) return;

    // Step 3: build new face list. For each merge component, walk the
    // boundary as an ordered HE loop and emit a single N-gon face.
    // Components of size 1 are kept unchanged.
    std::vector<Face>     newFaces;
    std::vector<HalfEdge> newHalfEdges;
    newFaces.reserve(compMembers.size());
    newHalfEdges.reserve(halfEdges_.size());

    // Map: old halfEdge index → new halfEdge index (NONE if discarded).
    std::vector<int32_t> heRemap(halfEdges_.size(), NONE);
    std::vector<int32_t> faceRemap(faces_.size(),   NONE);

    for (auto& kv : compMembers) {
        const std::vector<int32_t>& members = kv.second;

        if (members.size() == 1) {
            // Carry the face over verbatim.
            int32_t oldFi = members[0];
            int32_t newFi = (int32_t)newFaces.size();
            faceRemap[oldFi] = newFi;
            Face nf;
            nf.group = faces_[oldFi].group;
            int32_t start = faces_[oldFi].halfEdge;
            int32_t firstNew = NONE, prevNew = NONE;
            int32_t cur = start;
            do {
                int32_t newHi = (int32_t)newHalfEdges.size();
                newHalfEdges.push_back(halfEdges_[cur]);
                newHalfEdges.back().face = newFi;
                newHalfEdges.back().twin = NONE;     // re-linked at the end
                newHalfEdges.back().next = NONE;
                heRemap[cur] = newHi;
                if (prevNew != NONE) newHalfEdges[prevNew].next = newHi;
                else firstNew = newHi;
                prevNew = newHi;
                cur = halfEdges_[cur].next;
            } while (cur != start);
            newHalfEdges[prevNew].next = firstNew;
            nf.halfEdge = firstNew;
            newFaces.push_back(nf);
            continue;
        }

        // Multi-face component: walk the boundary.
        // Boundary HE = HE in component whose twin is NONE or in another component.
        std::unordered_set<int32_t> memberSet(members.begin(), members.end());
        std::vector<int32_t> boundary;
        for (int32_t fi : members) {
            int32_t start = faces_[fi].halfEdge;
            int32_t cur = start;
            do {
                int32_t tw = halfEdges_[cur].twin;
                bool isBd = (tw == NONE) ||
                            (memberSet.count(halfEdges_[tw].face) == 0);
                if (isBd) boundary.push_back(cur);
                cur = halfEdges_[cur].next;
            } while (cur != start);
        }
        if (boundary.empty()) continue;   // closed component → no boundary, skip

        // Order boundary into loop(s) using twin.next walks around the
        // shared destination vertex (same algorithm as findGroupBoundary).
        std::unordered_set<int32_t> bdSet(boundary.begin(), boundary.end());
        std::unordered_set<int32_t> visited;
        auto nextBd = [&](int32_t hi) -> int32_t {
            int32_t cur = halfEdges_[hi].next;
            size_t guard = 0;
            while (cur != NONE) {
                if (bdSet.count(cur)) return cur;
                int32_t tw = halfEdges_[cur].twin;
                if (tw == NONE) return NONE;
                cur = halfEdges_[tw].next;
                if (++guard > halfEdges_.size()) return NONE;
            }
            return NONE;
        };
        std::vector<std::vector<int32_t>> loops;
        for (int32_t start : boundary) {
            if (visited.count(start)) continue;
            std::vector<int32_t> loop;
            int32_t cur = start;
            while (cur != NONE && !visited.count(cur)) {
                visited.insert(cur);
                loop.push_back(cur);
                cur = nextBd(cur);
            }
            if (!loop.empty()) loops.push_back(std::move(loop));
        }

        // For each loop, emit one N-gon face. Most components have one loop
        // (= outer boundary). Multi-loop = inner holes which we encode as
        // separate faces with their own group tag for now.
        int32_t group = faces_[members[0]].group;
        for (const auto& loop : loops) {
            int32_t newFi = (int32_t)newFaces.size();
            int32_t firstNew = (int32_t)newHalfEdges.size();
            for (size_t i = 0; i < loop.size(); ++i) {
                HalfEdge nh;
                nh.origin = halfEdges_[loop[i]].origin;
                nh.face   = newFi;
                nh.twin   = NONE;
                nh.next   = (int32_t)(firstNew + ((i + 1) % loop.size()));
                newHalfEdges.push_back(nh);
                heRemap[loop[i]] = (int32_t)(firstNew + i);
            }
            Face nf;
            nf.halfEdge = firstNew;
            nf.group    = group;
            newFaces.push_back(nf);
            for (int32_t fi : members) faceRemap[fi] = newFi;
        }
    }

    // Step 4: re-link twins for HEs that survived. Two surviving HEs are
    // twins iff their original twins survived AND point to each other.
    for (size_t i = 0; i < halfEdges_.size(); ++i) {
        int32_t newHi = heRemap[i];
        if (newHi == NONE) continue;
        int32_t oldTw = halfEdges_[i].twin;
        if (oldTw == NONE) continue;
        int32_t newTw = heRemap[oldTw];
        if (newTw == NONE) continue;
        newHalfEdges[newHi].twin = newTw;
    }

    // Step 5: rebuild vertex outgoing-HE pointers — old indices are stale.
    for (Vertex& v : vertices_) v.halfEdge = NONE;
    for (size_t i = 0; i < newHalfEdges.size(); ++i) {
        int32_t oi = newHalfEdges[i].origin;
        if (oi >= 0 && oi < (int32_t)vertices_.size() &&
            vertices_[oi].halfEdge == NONE) {
            vertices_[oi].halfEdge = (int32_t)i;
        }
    }

    halfEdges_ = std::move(newHalfEdges);
    faces_     = std::move(newFaces);
}

// ===========================================================================
// Compact — drop unreferenced vertices and renumber
// ===========================================================================

void PolyMesh::compact() {
    std::vector<uint8_t> used(vertices_.size(), 0);
    for (const HalfEdge& he : halfEdges_) {
        if (he.origin >= 0 && he.origin < (int32_t)used.size()) used[he.origin] = 1;
    }
    std::vector<int32_t> remap(vertices_.size(), NONE);
    std::vector<Vertex> newVerts;
    newVerts.reserve(vertices_.size());
    for (size_t i = 0; i < vertices_.size(); ++i) {
        if (!used[i]) continue;
        remap[i] = (int32_t)newVerts.size();
        newVerts.push_back(vertices_[i]);
    }
    // Rewrite he origins.
    for (HalfEdge& he : halfEdges_) {
        if (he.origin >= 0 && he.origin < (int32_t)remap.size())
            he.origin = remap[he.origin];
    }
    // Reseed vertex halfEdge pointers.
    for (Vertex& v : newVerts) v.halfEdge = NONE;
    for (size_t i = 0; i < halfEdges_.size(); ++i) {
        int32_t oi = halfEdges_[i].origin;
        if (oi >= 0 && oi < (int32_t)newVerts.size() && newVerts[oi].halfEdge == NONE) {
            newVerts[oi].halfEdge = (int32_t)i;
        }
    }
    vertices_ = std::move(newVerts);
}

} // namespace bromesh
