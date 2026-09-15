#include "bromesh/manipulation/poly_mesh.h"
#include <bromath/vec.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace bromesh {

namespace {

constexpr float POS_QUANT = 1.0e5f;

struct PosKey {
    int64_t x, y, z;
    bool operator==(const PosKey& o) const { return x == o.x && y == o.y && z == o.z; }
};

struct PosKeyHash {
    size_t operator()(const PosKey& k) const {
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

} // namespace

// ===========================================================================
// Twin matching
// ===========================================================================

void PolyMesh::rematchTwins() {
    // Pass 1: vertex-identity pairing. Key = (origin, destination)
    std::unordered_map<uint64_t, int32_t> byIdDir;
    byIdDir.reserve(halfEdges_.size());
    auto idKey = [](int32_t a, int32_t b) -> uint64_t {
        return ((uint64_t)(uint32_t)a << 32) | (uint64_t)(uint32_t)b;
    };

    for (size_t i = 0; i < halfEdges_.size(); ++i) {
        if (!isLiveHalfEdge((int32_t)i)) continue;
        int32_t orig = halfEdges_[i].origin;
        int32_t dest = halfEdges_[halfEdges_[i].next].origin;
        if (orig == NONE || dest == NONE) continue;
        byIdDir.emplace(idKey(orig, dest), (int32_t)i);
    }

    for (size_t i = 0; i < halfEdges_.size(); ++i) {
        if (!isLiveHalfEdge((int32_t)i)) continue;
        if (halfEdges_[i].twin != NONE) continue;
        int32_t orig = halfEdges_[i].origin;
        int32_t dest = halfEdges_[halfEdges_[i].next].origin;
        if (orig == NONE || dest == NONE) continue;
        auto it = byIdDir.find(idKey(dest, orig));
        if (it != byIdDir.end()) {
            int32_t tj = it->second;
            if (halfEdges_[tj].twin == NONE) {
                halfEdges_[i].twin = tj;
                halfEdges_[tj].twin = (int32_t)i;
            }
        }
    }

    // Pass 2: position-based pairing for hard-edge seams
    struct DirectedEdgePos {
        PosKey a, b;
        bool operator==(const DirectedEdgePos& o) const { return a == o.a && b == o.b; }
    };
    struct DirectedEdgePosHash {
        size_t operator()(const DirectedEdgePos& e) const {
            PosKeyHash h;
            return h(e.a) ^ (h(e.b) + 0x9e3779b9 + (h(e.a) << 6) + (h(e.a) >> 2));
        }
    };

    auto posKeyOfVert = [&](int32_t vi) -> PosKey {
        const Vertex& v = vertices_[vi];
        return keyOf(v.x, v.y, v.z);
    };

    std::unordered_map<DirectedEdgePos, int32_t, DirectedEdgePosHash> byPosDir;
    for (size_t i = 0; i < halfEdges_.size(); ++i) {
        if (!isLiveHalfEdge((int32_t)i)) continue;
        if (halfEdges_[i].twin != NONE) continue;
        int32_t orig = halfEdges_[i].origin;
        int32_t dest = halfEdges_[halfEdges_[i].next].origin;
        if (orig == NONE || dest == NONE) continue;
        byPosDir.emplace(DirectedEdgePos{posKeyOfVert(orig), posKeyOfVert(dest)}, (int32_t)i);
    }

    for (size_t i = 0; i < halfEdges_.size(); ++i) {
        if (!isLiveHalfEdge((int32_t)i)) continue;
        if (halfEdges_[i].twin != NONE) continue;
        int32_t orig = halfEdges_[i].origin;
        int32_t dest = halfEdges_[halfEdges_[i].next].origin;
        if (orig == NONE || dest == NONE) continue;
        auto it = byPosDir.find(DirectedEdgePos{posKeyOfVert(dest), posKeyOfVert(orig)});
        if (it == byPosDir.end()) continue;
        int32_t tj = it->second;
        if (halfEdges_[tj].twin != NONE) continue;
        halfEdges_[i].twin = tj;
        halfEdges_[tj].twin = (int32_t)i;
    }
}

// ===========================================================================
// Surgery — translate / extrude / inset
// ===========================================================================

void PolyMesh::translateVertex(int32_t vi, const float offset[3]) {
    if (vi < 0 || vi >= (int32_t)vertices_.size()) return;
    vertices_[vi].x += offset[0];
    vertices_[vi].y += offset[1];
    vertices_[vi].z += offset[2];
}

void PolyMesh::translateFace(int32_t faceIdx, const float offset[3]) {
    auto verts = faceVertices(faceIdx);
    std::unordered_set<int32_t> seen;
    for (int32_t vi : verts) {
        if (seen.insert(vi).second) translateVertex(vi, offset);
    }
}

void PolyMesh::translateFaceWithRing(int32_t faceIdx, const float offset[3]) {
    auto verts = faceVertices(faceIdx);
    if (verts.empty()) return;

    std::unordered_set<int32_t> allVerts(verts.begin(), verts.end());
    for (int32_t seedVi : verts) {
        int32_t startH = vertices_[seedVi].halfEdge;
        if (startH == NONE) continue;
        int32_t curH = startH;
        size_t guard = 0;
        do {
            allVerts.insert(halfEdges_[curH].origin);
            int32_t tw = halfEdges_[curH].twin;
            if (tw == NONE) break;
            curH = halfEdges_[tw].next;
            if (++guard > halfEdges_.size()) break;
        } while (curH != startH && curH != NONE);
    }

    std::vector<PosKey> targets;
    targets.reserve(verts.size());
    for (int32_t vi : verts) {
        const auto& v = vertices_[vi];
        targets.push_back(keyOf(v.x, v.y, v.z));
    }
    for (size_t vi = 0; vi < vertices_.size(); ++vi) {
        const auto& v = vertices_[vi];
        PosKey k = keyOf(v.x, v.y, v.z);
        for (const auto& tk : targets) {
            if (k == tk) { allVerts.insert((int32_t)vi); break; }
        }
    }

    for (int32_t vi : allVerts) translateVertex(vi, offset);
}

PolyMesh::ExtrudeResult PolyMesh::extrudeFace(int32_t faceIdx,
                                              const float offset[3],
                                              bool withBackFace,
                                              int32_t bridgeGroup,
                                              int32_t backGroup) {
    ExtrudeResult out;
    if (faceIdx < 0 || faceIdx >= (int32_t)faces_.size()) return out;

    std::vector<int32_t> bdHEs = faceHalfEdges(faceIdx);
    if (bdHEs.size() < 3) return out;

    struct BoundaryRec {
        int32_t he, oldA, oldB, newA = NONE, newB = NONE, adjGroup = -1;
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

    for (int32_t hi : bdHEs) {
        halfEdges_[hi].origin = dupMap.at(halfEdges_[hi].origin);
        if (vertices_[halfEdges_[hi].origin].halfEdge == NONE) {
            vertices_[halfEdges_[hi].origin].halfEdge = hi;
        }
    }

    for (int32_t hi : bdHEs) {
        int32_t tw = halfEdges_[hi].twin;
        if (tw != NONE) {
            halfEdges_[tw].twin = NONE;
            halfEdges_[hi].twin = NONE;
        }
    }

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
        std::vector<int32_t> back;
        back.reserve(recs.size());
        for (auto it = recs.rbegin(); it != recs.rend(); ++it) back.push_back(it->oldA);
        out.backFace = addFace(back, backGroup);
    }

    rematchTwins();
    return out;
}

PolyMesh::InsetResult PolyMesh::insetFace(int32_t faceIdx, float amount, bool asRatio, int32_t bridgeGroup) {
    InsetResult out;
    if (!isLiveFace(faceIdx) || amount <= 0.0f) return out;
    if (asRatio && amount >= 1.0f) return out;

    const std::vector<int32_t> outerVerts = faceVertices(faceIdx);
    const std::vector<int32_t> bdHEs = faceHalfEdges(faceIdx);
    const size_t N = outerVerts.size();
    if (N < 3 || bdHEs.size() != N) return out;

    float normArr[3];
    computeFaceNormal(faceIdx, normArr);
    bromath::Vec3 N_norm(normArr[0], normArr[1], normArr[2]);
    if (bromath::vlen2(N_norm) < 1e-12f) return out;

    bromath::Vec3 centroid(0.0f, 0.0f, 0.0f);
    for (int32_t vi : outerVerts) {
        centroid += bromath::Vec3(vertices_[vi].x, vertices_[vi].y, vertices_[vi].z);
    }
    centroid = centroid * (1.0f / static_cast<float>(N));

    std::vector<int32_t> innerVerts;
    innerVerts.reserve(N);
    for (size_t i = 0; i < N; ++i) {
        int32_t vi = outerVerts[i];
        bromath::Vec3 p(vertices_[vi].x, vertices_[vi].y, vertices_[vi].z);
        bromath::Vec3 pInner;
        if (asRatio) {
            pInner = p + (centroid - p) * amount;
        } else {
            int32_t prevVi = outerVerts[(i + N - 1) % N];
            int32_t nextVi = outerVerts[(i + 1) % N];
            bromath::Vec3 pPrev(vertices_[prevVi].x, vertices_[prevVi].y, vertices_[prevVi].z);
            bromath::Vec3 pNext(vertices_[nextVi].x, vertices_[nextVi].y, vertices_[nextVi].z);

            bromath::Vec3 t0 = bromath::vnorm(p - pPrev);
            bromath::Vec3 t1 = bromath::vnorm(pNext - p);
            bromath::Vec3 n0 = bromath::vcross(N_norm, t0);
            bromath::Vec3 n1 = bromath::vcross(N_norm, t1);
            bromath::Vec3 bisector = n0 + n1;
            float bLen = bromath::vlen(bisector);
            if (bLen > 1e-6f) {
                bisector = bisector / bLen;
                float denom = bromath::vdot(bisector, n0);
                float dist = (denom > 0.25f) ? std::min(amount / denom, amount * 4.0f) : (amount * 4.0f);
                float maxDist = bromath::vdist(p, centroid) * 0.99f;
                if (dist > maxDist && maxDist > 0.0f) dist = maxDist;
                pInner = p + bisector * dist;
            } else {
                bromath::Vec3 toC = centroid - p;
                float dC = bromath::vlen(toC);
                float dist = std::min(amount, dC * 0.99f);
                pInner = p + bromath::vnormOr(toC, bromath::Vec3(0, 0, 0)) * dist;
            }
        }
        int32_t inVi = addVertex(pInner.x, pInner.y, pInner.z);
        innerVerts.push_back(inVi);
    }

    for (int32_t hi : bdHEs) {
        int32_t tw = halfEdges_[hi].twin;
        if (tw != NONE) {
            halfEdges_[tw].twin = NONE;
            halfEdges_[hi].twin = NONE;
        }
    }

    for (size_t i = 0; i < N; ++i) {
        int32_t hi = bdHEs[i];
        halfEdges_[hi].origin = innerVerts[i];
        vertices_[innerVerts[i]].halfEdge = hi;
    }

    for (int32_t vi : outerVerts) {
        if (vertices_[vi].halfEdge != NONE &&
            halfEdges_[vertices_[vi].halfEdge].origin != vi) {
            vertices_[vi].halfEdge = NONE;
        }
    }

    std::vector<int32_t> bridgeFaces;
    bridgeFaces.reserve(N);
    for (size_t i = 0; i < N; ++i) {
        int32_t o0 = outerVerts[i];
        int32_t o1 = outerVerts[(i + 1) % N];
        int32_t i1 = innerVerts[(i + 1) % N];
        int32_t i0 = innerVerts[i];
        int32_t bFi = addFace({o0, o1, i1, i0}, bridgeGroup);
        bridgeFaces.push_back(bFi);
    }

    for (int32_t vi : outerVerts) {
        if (vertices_[vi].halfEdge == NONE) {
            for (size_t h = 0; h < halfEdges_.size(); ++h) {
                if (isLiveHalfEdge((int32_t)h) && halfEdges_[h].origin == vi) {
                    vertices_[vi].halfEdge = (int32_t)h;
                    break;
                }
            }
        }
    }

    rematchTwins();

    out.innerFace = faceIdx;
    out.innerVerts = std::move(innerVerts);
    out.bridgeFaces = std::move(bridgeFaces);
    return out;
}

// ===========================================================================
// Surgery — split / flip / collapse / delete
// ===========================================================================

int32_t PolyMesh::splitEdge(int32_t hi, const float* posOptional) {
    if (!isLiveHalfEdge(hi)) return NONE;

    const int32_t h1 = hi, h2 = halfEdges_[h1].next, h3 = halfEdges_[h2].next;
    if (halfEdges_[h3].next != h1) return NONE;
    const int32_t a = halfEdges_[h1].origin, b = halfEdges_[h2].origin, c = halfEdges_[h3].origin;
    const int32_t F1 = halfEdges_[h1].face;

    const int32_t t1 = halfEdges_[h1].twin;
    const bool boundary = (t1 == NONE);
    int32_t t2 = NONE, t3 = NONE, d = NONE, F2 = NONE;
    if (!boundary) {
        t2 = halfEdges_[t1].next;
        t3 = halfEdges_[t2].next;
        if (halfEdges_[t3].next != t1) return NONE;
        d = halfEdges_[t3].origin;
        F2 = halfEdges_[t1].face;
    }

    float mx, my, mz;
    if (posOptional) {
        mx = posOptional[0]; my = posOptional[1]; mz = posOptional[2];
    } else {
        const Vertex& va = vertices_[a];
        const Vertex& vb = vertices_[b];
        mx = (va.x + vb.x) * 0.5f; my = (va.y + vb.y) * 0.5f; mz = (va.z + vb.z) * 0.5f;
    }
    const int32_t M = addVertex(mx, my, mz);

    const int32_t h_MB = (int32_t)halfEdges_.size(); halfEdges_.push_back({});
    const int32_t h_AM_to_c = (int32_t)halfEdges_.size(); halfEdges_.push_back({});
    const int32_t h_C_to_M = (int32_t)halfEdges_.size(); halfEdges_.push_back({});
    const int32_t F1_new = (int32_t)faces_.size();
    faces_.push_back(Face{ h_MB, faces_[F1].group });

    halfEdges_[h1].origin = a; halfEdges_[h1].next = h_AM_to_c; halfEdges_[h1].face = F1;
    halfEdges_[h_AM_to_c] = HalfEdge{ M, NONE, h3, F1 };
    faces_[F1].halfEdge = h1;

    halfEdges_[h_MB] = HalfEdge{ M, NONE, h2, F1_new };
    halfEdges_[h2].next = h_C_to_M; halfEdges_[h2].face = F1_new;
    halfEdges_[h_C_to_M] = HalfEdge{ c, NONE, h_MB, F1_new };

    halfEdges_[h_AM_to_c].twin = h_C_to_M;
    halfEdges_[h_C_to_M].twin = h_AM_to_c;

    if (!boundary) {
        const int32_t h_M_to_d = (int32_t)halfEdges_.size(); halfEdges_.push_back({});
        const int32_t h_M_to_a = (int32_t)halfEdges_.size(); halfEdges_.push_back({});
        const int32_t h_D_to_M = (int32_t)halfEdges_.size(); halfEdges_.push_back({});
        const int32_t F2_new = (int32_t)faces_.size();
        faces_.push_back(Face{ h_M_to_a, faces_[F2].group });

        halfEdges_[t1].origin = b; halfEdges_[t1].next = h_M_to_d; halfEdges_[t1].face = F2;
        halfEdges_[h_M_to_d] = HalfEdge{ M, NONE, t3, F2 };
        faces_[F2].halfEdge = t1;

        halfEdges_[h_M_to_a] = HalfEdge{ M, NONE, t2, F2_new };
        halfEdges_[t2].next = h_D_to_M; halfEdges_[t2].face = F2_new;
        halfEdges_[h_D_to_M] = HalfEdge{ d, NONE, h_M_to_a, F2_new };

        halfEdges_[h_M_to_d].twin = h_D_to_M;
        halfEdges_[h_D_to_M].twin = h_M_to_d;

        halfEdges_[h1].twin = h_M_to_a;
        halfEdges_[h_M_to_a].twin = h1;
        halfEdges_[t1].twin = h_MB;
        halfEdges_[h_MB].twin = t1;
    } else {
        halfEdges_[h1].twin = NONE;
        halfEdges_[h_MB].twin = NONE;
    }

    vertices_[a].halfEdge = h1; vertices_[b].halfEdge = h2;
    vertices_[c].halfEdge = h3;
    if (!boundary) vertices_[d].halfEdge = t3;
    vertices_[M].halfEdge = h_AM_to_c;

    return M;
}

bool PolyMesh::flipEdge(int32_t hi) {
    if (!isLiveHalfEdge(hi)) return false;
    const int32_t h1 = hi, t1 = halfEdges_[h1].twin;
    if (t1 == NONE) return false;

    const int32_t h2 = halfEdges_[h1].next, h3 = halfEdges_[h2].next;
    if (halfEdges_[h3].next != h1) return false;

    const int32_t t2 = halfEdges_[t1].next, t3 = halfEdges_[t2].next;
    if (halfEdges_[t3].next != t1) return false;

    const int32_t a = halfEdges_[h1].origin, b = halfEdges_[h2].origin;
    const int32_t c = halfEdges_[h3].origin, d = halfEdges_[t3].origin;
    if (c == d) return false;

    const int32_t F1 = halfEdges_[h1].face, F2 = halfEdges_[t1].face;

    for (size_t i = 0; i < halfEdges_.size(); ++i) {
        if (!isLiveHalfEdge((int32_t)i)) continue;
        const HalfEdge& he = halfEdges_[i];
        if ((he.origin == c && halfEdges_[he.next].origin == d) ||
            (he.origin == d && halfEdges_[he.next].origin == c)) return false;
    }

    halfEdges_[h1].origin = d; halfEdges_[h1].next = h3; halfEdges_[h1].face = F1;
    halfEdges_[t1].origin = c; halfEdges_[t1].next = t3; halfEdges_[t1].face = F2;

    halfEdges_[t2].next = h1; halfEdges_[t2].face = F1; halfEdges_[h3].next = t2;
    halfEdges_[h2].next = t1; halfEdges_[h2].face = F2; halfEdges_[t3].next = h2;

    faces_[F1].halfEdge = t2; faces_[F2].halfEdge = h2;
    vertices_[a].halfEdge = t2; vertices_[b].halfEdge = h2;
    vertices_[c].halfEdge = h3; vertices_[d].halfEdge = t3;

    return true;
}

bool PolyMesh::collapseEdge(int32_t hi, const float* posOptional) {
    if (!isLiveHalfEdge(hi)) return false;

    const int32_t h1 = hi, t1 = halfEdges_[h1].twin;
    const bool boundary = (t1 == NONE);

    const int32_t h2 = halfEdges_[h1].next, h3 = halfEdges_[h2].next;
    if (halfEdges_[h3].next != h1) return false;
    const int32_t a = halfEdges_[h1].origin, b = halfEdges_[h2].origin, c = halfEdges_[h3].origin;
    const int32_t F1 = halfEdges_[h1].face;

    int32_t t2 = NONE, t3 = NONE, d = NONE, F2 = NONE;
    if (!boundary) {
        t2 = halfEdges_[t1].next;
        t3 = halfEdges_[t2].next;
        if (halfEdges_[t3].next != t1) return false;
        d = halfEdges_[t3].origin;
        F2 = halfEdges_[t1].face;
    }

    if (!boundary) {
        if (isBoundaryVertex(a) || isBoundaryVertex(b)) return false;
    }

    auto neighborSet = [&](int32_t v) {
        std::unordered_set<int32_t> s;
        for (size_t i = 0; i < halfEdges_.size(); ++i) {
            if (!isLiveHalfEdge((int32_t)i)) continue;
            if (halfEdges_[i].origin != v) continue;
            s.insert(halfEdges_[halfEdges_[i].next].origin);
        }
        return s;
    };
    auto nA = neighborSet(a), nB = neighborSet(b);
    for (int32_t n : nB) {
        if (n == a || n == b) continue;
        if (!nA.count(n)) continue;
        if (n == c) continue;
        if (!boundary && n == d) continue;
        return false;
    }

    auto stitch = [&](int32_t wingA, int32_t wingB) {
        int32_t extA = halfEdges_[wingA].twin;
        int32_t extB = halfEdges_[wingB].twin;
        if (extA != NONE) halfEdges_[extA].twin = extB;
        if (extB != NONE) halfEdges_[extB].twin = extA;
    };
    stitch(h2, h3);
    if (!boundary) stitch(t2, t3);

    if (posOptional) {
        vertices_[a].x = posOptional[0]; vertices_[a].y = posOptional[1]; vertices_[a].z = posOptional[2];
    } else {
        vertices_[a].x = 0.5f * (vertices_[a].x + vertices_[b].x);
        vertices_[a].y = 0.5f * (vertices_[a].y + vertices_[b].y);
        vertices_[a].z = 0.5f * (vertices_[a].z + vertices_[b].z);
    }

    for (size_t i = 0; i < halfEdges_.size(); ++i) {
        if (!isLiveHalfEdge((int32_t)i)) continue;
        if (halfEdges_[i].face == F1) continue;
        if (!boundary && halfEdges_[i].face == F2) continue;
        if (halfEdges_[i].origin == b) halfEdges_[i].origin = a;
    }

    auto tombstoneFace = [&](int32_t fi, std::initializer_list<int32_t> ring) {
        for (int32_t hh : ring) {
            halfEdges_[hh].origin = NONE; halfEdges_[hh].next = NONE;
            halfEdges_[hh].face = NONE;   halfEdges_[hh].twin = NONE;
        }
        faces_[fi].halfEdge = NONE;
    };
    tombstoneFace(F1, { h1, h2, h3 });
    if (!boundary) tombstoneFace(F2, { t1, t2, t3 });

    vertices_[a].halfEdge = NONE;
    for (size_t i = 0; i < halfEdges_.size(); ++i) {
        if (!isLiveHalfEdge((int32_t)i)) continue;
        if (halfEdges_[i].origin == a) { vertices_[a].halfEdge = (int32_t)i; break; }
    }
    vertices_[b].halfEdge = NONE;

    auto reseedVertex = [&](int32_t vi) {
        if (vi < 0 || vi >= (int32_t)vertices_.size()) return;
        if (isLiveHalfEdge(vertices_[vi].halfEdge) &&
            halfEdges_[vertices_[vi].halfEdge].origin == vi) return;
        vertices_[vi].halfEdge = NONE;
        for (size_t i = 0; i < halfEdges_.size(); ++i) {
            if (!isLiveHalfEdge((int32_t)i)) continue;
            if (halfEdges_[i].origin == vi) { vertices_[vi].halfEdge = (int32_t)i; break; }
        }
    };
    reseedVertex(c);
    if (!boundary) reseedVertex(d);

    return true;
}

void PolyMesh::deleteFace(int32_t fi) {
    if (!isLiveFace(fi)) return;

    std::vector<int32_t> ring;
    ring.reserve(8);
    int32_t start = faces_[fi].halfEdge, cur = start;
    size_t guard = 0;
    do {
        ring.push_back(cur);
        cur = halfEdges_[cur].next;
        if (++guard > halfEdges_.size()) break;
    } while (cur != start && cur != NONE);

    for (int32_t hi : ring) {
        int32_t tw = halfEdges_[hi].twin;
        if (tw != NONE && tw >= 0 && tw < (int32_t)halfEdges_.size()) {
            halfEdges_[tw].twin = NONE;
        }
    }

    std::vector<int32_t> orphanedOrigins;
    orphanedOrigins.reserve(ring.size());
    for (int32_t hi : ring) {
        orphanedOrigins.push_back(halfEdges_[hi].origin);
        halfEdges_[hi].origin = NONE; halfEdges_[hi].next = NONE;
        halfEdges_[hi].face = NONE;   halfEdges_[hi].twin = NONE;
    }
    faces_[fi].halfEdge = NONE;

    for (int32_t vi : orphanedOrigins) {
        if (vi < 0 || vi >= (int32_t)vertices_.size()) continue;
        if (isLiveHalfEdge(vertices_[vi].halfEdge)) continue;
        vertices_[vi].halfEdge = NONE;
        for (size_t i = 0; i < halfEdges_.size(); ++i) {
            if (halfEdges_[i].face == NONE) continue;
            if (halfEdges_[i].origin == vi) {
                vertices_[vi].halfEdge = (int32_t)i;
                break;
            }
        }
    }
}

// ===========================================================================
// mergeFacesByGroup — collapse coplanar tri groups into N-gon faces
// ===========================================================================

void PolyMesh::mergeFacesByGroup() {
    if (faces_.empty()) return;

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
        int32_t fa = he.face, fb = halfEdges_[he.twin].face;
        if (fa < 0 || fb < 0) continue;
        if (faces_[fa].group != faces_[fb].group) continue;
        if (faces_[fa].group < 0) continue;
        unite(fa, fb);
    }

    std::unordered_map<int32_t, std::vector<int32_t>> compMembers;
    for (size_t i = 0; i < faces_.size(); ++i) {
        compMembers[find((int32_t)i)].push_back((int32_t)i);
    }
    bool anyMerge = false;
    for (auto& kv : compMembers) if (kv.second.size() > 1) { anyMerge = true; break; }
    if (!anyMerge) return;

    std::vector<Face> newFaces;
    std::vector<HalfEdge> newHalfEdges;
    newFaces.reserve(compMembers.size());
    newHalfEdges.reserve(halfEdges_.size());

    std::vector<int32_t> heRemap(halfEdges_.size(), NONE);
    std::vector<int32_t> faceRemap(faces_.size(), NONE);

    for (auto& kv : compMembers) {
        const std::vector<int32_t>& members = kv.second;

        if (members.size() == 1) {
            int32_t oldFi = members[0], newFi = (int32_t)newFaces.size();
            faceRemap[oldFi] = newFi;
            Face nf;
            nf.group = faces_[oldFi].group;
            int32_t start = faces_[oldFi].halfEdge, firstNew = NONE, prevNew = NONE, cur = start;
            do {
                int32_t newHi = (int32_t)newHalfEdges.size();
                newHalfEdges.push_back(halfEdges_[cur]);
                newHalfEdges.back().face = newFi;
                newHalfEdges.back().twin = NONE;
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

        std::unordered_set<int32_t> memberSet(members.begin(), members.end());
        std::vector<int32_t> boundary;
        for (int32_t fi : members) {
            int32_t start = faces_[fi].halfEdge, cur = start;
            do {
                int32_t tw = halfEdges_[cur].twin;
                bool isBd = (tw == NONE) || (memberSet.count(halfEdges_[tw].face) == 0);
                if (isBd) boundary.push_back(cur);
                cur = halfEdges_[cur].next;
            } while (cur != start);
        }
        if (boundary.empty()) continue;

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

        int32_t grp = faces_[members[0]].group;
        for (const auto& loop : loops) {
            if (loop.size() < 3) continue;
            int32_t newFi = (int32_t)newFaces.size();
            for (int32_t oldFi : members) faceRemap[oldFi] = newFi;
            Face nf;
            nf.group = grp;
            int32_t firstNew = NONE, prevNew = NONE;
            for (int32_t oldHi : loop) {
                int32_t newHi = (int32_t)newHalfEdges.size();
                HalfEdge he;
                he.origin = halfEdges_[oldHi].origin;
                he.face = newFi;
                he.twin = NONE;
                he.next = NONE;
                newHalfEdges.push_back(he);
                heRemap[oldHi] = newHi;
                if (prevNew != NONE) newHalfEdges[prevNew].next = newHi;
                else firstNew = newHi;
                prevNew = newHi;
            }
            newHalfEdges[prevNew].next = firstNew;
            nf.halfEdge = firstNew;
            newFaces.push_back(nf);
        }
    }

    for (size_t oldHi = 0; oldHi < halfEdges_.size(); ++oldHi) {
        int32_t newHi = heRemap[oldHi];
        if (newHi == NONE) continue;
        int32_t oldTw = halfEdges_[oldHi].twin;
        if (oldTw != NONE) {
            int32_t newTw = heRemap[oldTw];
            if (newTw != NONE) newHalfEdges[newHi].twin = newTw;
        }
    }

    for (auto& v : vertices_) v.halfEdge = NONE;
    for (size_t i = 0; i < newHalfEdges.size(); ++i) {
        int32_t oi = newHalfEdges[i].origin;
        if (oi >= 0 && oi < (int32_t)vertices_.size() && vertices_[oi].halfEdge == NONE) {
            vertices_[oi].halfEdge = (int32_t)i;
        }
    }

    halfEdges_ = std::move(newHalfEdges);
    faces_ = std::move(newFaces);
}

} // namespace bromesh
