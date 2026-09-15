#include "bromesh/manipulation/poly_mesh.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

namespace bromesh {

namespace {

static void normalize3(float v[3]) {
    float L = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if (L > 1.0e-20f) { v[0] /= L; v[1] /= L; v[2] /= L; }
    else              { v[0] = v[1] = v[2] = 0.0f; }
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
        if (count > (int)halfEdges_.size()) return 0;
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

bool PolyMesh::isBoundaryVertex(int32_t vi) const {
    if (vi < 0 || vi >= (int32_t)vertices_.size()) return false;
    for (size_t i = 0; i < halfEdges_.size(); ++i) {
        const HalfEdge& he = halfEdges_[i];
        if (he.face == NONE) continue;
        if (he.twin != NONE) continue;
        if (he.origin == vi) return true;
        if (halfEdges_[he.next].origin == vi) return true;
    }
    return false;
}

// ===========================================================================
// Boundary discovery
// ===========================================================================

std::vector<std::vector<int32_t>> PolyMesh::findFaceBoundary(int faceIdx) const {
    std::vector<std::vector<int32_t>> loops;
    if (faceIdx < 0 || faceIdx >= (int)faces_.size()) return loops;
    std::vector<int32_t> loop;
    int32_t start = faces_[faceIdx].halfEdge, cur = start;
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

    std::unordered_set<int32_t> visited;
    for (int32_t start : seeds) {
        if (visited.count(start)) continue;
        std::vector<int32_t> loop;
        int32_t cur = start;
        while (cur != NONE && !visited.count(cur)) {
            visited.insert(cur);
            loop.push_back(cur);
            int32_t next = halfEdges_[cur].next, chosen = NONE;
            size_t guard = 0;
            while (next != NONE) {
                if (faces_[halfEdges_[next].face].group == groupId && isGroupBoundary(next)) {
                    chosen = next;
                    break;
                }
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
        int32_t h1 = h0 + 1, h2 = h0 + 2;
        pm.halfEdges_.push_back(HalfEdge{(int32_t)indices[t*3 + 0], NONE, h1, fi});
        pm.halfEdges_.push_back(HalfEdge{(int32_t)indices[t*3 + 1], NONE, h2, fi});
        pm.halfEdges_.push_back(HalfEdge{(int32_t)indices[t*3 + 2], NONE, h0, fi});
        f.halfEdge = h0;
        pm.faces_.push_back(f);
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
        const uint32_t end = polyOffsets[f + 1];
        if (end - start < 3) continue;
        std::vector<int32_t> verts;
        verts.reserve(end - start);
        for (uint32_t k = start; k < end; ++k) verts.push_back((int32_t)polyVerts[k]);
        pm.addFace(verts, f < faceGroups.size() ? faceGroups[f] : -1);
    }
    pm.rematchTwins();
    return pm;
}

int32_t PolyMesh::addVertex(float x, float y, float z) {
    vertices_.push_back(Vertex{x, y, z, NONE});
    return (int32_t)vertices_.size() - 1;
}

int32_t PolyMesh::addFace(const std::vector<int32_t>& verts, int32_t group) {
    if (verts.size() < 3) return NONE;
    const int32_t fi = (int32_t)faces_.size();
    const int32_t hBase = (int32_t)halfEdges_.size();
    const int32_t N = (int32_t)verts.size();
    for (int32_t k = 0; k < N; ++k) {
        HalfEdge he;
        he.origin = verts[k];
        he.twin = NONE;
        he.next = hBase + ((k + 1) % N);
        he.face = fi;
        halfEdges_.push_back(he);
        if (vertices_[verts[k]].halfEdge == NONE) {
            vertices_[verts[k]].halfEdge = hBase + k;
        }
    }
    Face f;
    f.halfEdge = hBase;
    f.group = group;
    faces_.push_back(f);
    return fi;
}

// ===========================================================================
// Validation
// ===========================================================================

PolyMesh::Validation PolyMesh::validate() const {
    Validation r;
    for (size_t i = 0; i < halfEdges_.size(); ++i) {
        const HalfEdge& he = halfEdges_[i];
        if (he.face == NONE) continue;
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
        if (f.halfEdge == NONE) continue;
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
// Compact — drop unreferenced vertices, dead faces, dead half-edges
// ===========================================================================

void PolyMesh::compact() {
    std::vector<int32_t> heRemap(halfEdges_.size(), NONE);
    std::vector<HalfEdge> newHe;
    newHe.reserve(halfEdges_.size());
    for (size_t i = 0; i < halfEdges_.size(); ++i) {
        if (halfEdges_[i].face == NONE) continue;
        heRemap[i] = (int32_t)newHe.size();
        newHe.push_back(halfEdges_[i]);
    }

    std::vector<int32_t> faceRemap(faces_.size(), NONE);
    std::vector<Face> newFaces;
    newFaces.reserve(faces_.size());
    for (size_t i = 0; i < faces_.size(); ++i) {
        if (faces_[i].halfEdge == NONE) continue;
        faceRemap[i] = (int32_t)newFaces.size();
        newFaces.push_back(faces_[i]);
    }

    for (HalfEdge& he : newHe) {
        if (he.next >= 0 && he.next < (int32_t)heRemap.size())
            he.next = heRemap[he.next];
        if (he.twin >= 0 && he.twin < (int32_t)heRemap.size())
            he.twin = heRemap[he.twin];
        if (he.face >= 0 && he.face < (int32_t)faceRemap.size())
            he.face = faceRemap[he.face];
    }
    for (Face& f : newFaces) {
        if (f.halfEdge >= 0 && f.halfEdge < (int32_t)heRemap.size())
            f.halfEdge = heRemap[f.halfEdge];
    }

    halfEdges_ = std::move(newHe);
    faces_     = std::move(newFaces);

    std::vector<uint8_t> used(vertices_.size(), 0);
    for (const HalfEdge& he : halfEdges_) {
        if (he.origin >= 0 && he.origin < (int32_t)used.size()) used[he.origin] = 1;
    }
    std::vector<int32_t> vremap(vertices_.size(), NONE);
    std::vector<Vertex> newVerts;
    newVerts.reserve(vertices_.size());
    for (size_t i = 0; i < vertices_.size(); ++i) {
        if (!used[i]) continue;
        vremap[i] = (int32_t)newVerts.size();
        newVerts.push_back(vertices_[i]);
    }
    for (HalfEdge& he : halfEdges_) {
        if (he.origin >= 0 && he.origin < (int32_t)vremap.size())
            he.origin = vremap[he.origin];
    }
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
