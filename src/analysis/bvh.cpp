#include "bromesh/analysis/bvh.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>

namespace bromesh {

// --- Local helpers -----------------------------------------------------------

namespace {

inline void triBounds(const MeshData& mesh, uint32_t tri,
                      float outMin[3], float outMax[3],
                      float outCentroid[3])
{
    const uint32_t i0 = mesh.indices[tri * 3 + 0];
    const uint32_t i1 = mesh.indices[tri * 3 + 1];
    const uint32_t i2 = mesh.indices[tri * 3 + 2];
    const float* p0 = &mesh.positions[i0 * 3];
    const float* p1 = &mesh.positions[i1 * 3];
    const float* p2 = &mesh.positions[i2 * 3];
    for (int a = 0; a < 3; ++a) {
        float lo = std::min(p0[a], std::min(p1[a], p2[a]));
        float hi = std::max(p0[a], std::max(p1[a], p2[a]));
        outMin[a] = lo;
        outMax[a] = hi;
        outCentroid[a] = (p0[a] + p1[a] + p2[a]) * (1.0f / 3.0f);
    }
}

inline void expandBBox(float dstMin[3], float dstMax[3],
                       const float srcMin[3], const float srcMax[3])
{
    for (int a = 0; a < 3; ++a) {
        if (srcMin[a] < dstMin[a]) dstMin[a] = srcMin[a];
        if (srcMax[a] > dstMax[a]) dstMax[a] = srcMax[a];
    }
}

inline void expandPoint(float dstMin[3], float dstMax[3], const float p[3])
{
    for (int a = 0; a < 3; ++a) {
        if (p[a] < dstMin[a]) dstMin[a] = p[a];
        if (p[a] > dstMax[a]) dstMax[a] = p[a];
    }
}

// Möller-Trumbore — returns true with (t, u, v) if the ray hits the triangle
// at positive distance. u,v are the barycentric coords for v1,v2 (v0 weight
// implicit = 1 - u - v).
inline bool rayTriangle(
    const float* orig, const float* dir,
    const float* v0, const float* v1, const float* v2,
    float& t, float& u, float& v)
{
    float e1[3] = { v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2] };
    float e2[3] = { v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2] };

    float h[3] = {
        dir[1] * e2[2] - dir[2] * e2[1],
        dir[2] * e2[0] - dir[0] * e2[2],
        dir[0] * e2[1] - dir[1] * e2[0]
    };
    float a = e1[0] * h[0] + e1[1] * h[1] + e1[2] * h[2];
    if (std::fabs(a) < 1e-10f) return false;

    float f = 1.0f / a;
    float s[3] = { orig[0] - v0[0], orig[1] - v0[1], orig[2] - v0[2] };
    u = f * (s[0] * h[0] + s[1] * h[1] + s[2] * h[2]);
    if (u < 0.0f || u > 1.0f) return false;

    float q[3] = {
        s[1] * e1[2] - s[2] * e1[1],
        s[2] * e1[0] - s[0] * e1[2],
        s[0] * e1[1] - s[1] * e1[0]
    };
    v = f * (dir[0] * q[0] + dir[1] * q[1] + dir[2] * q[2]);
    if (v < 0.0f || u + v > 1.0f) return false;

    t = f * (e2[0] * q[0] + e2[1] * q[1] + e2[2] * q[2]);
    return t > 0.0f;
}

// Kay-Kajiya slab test. invDir is 1/direction (component-wise). Returns the
// entry distance tMin if the ray hits the box within [0, tMax], else a
// negative sentinel. Handles infinities (axis-aligned rays) correctly because
// IEEE 0*inf produces NaN which the min/max branches tolerate: we use a form
// that avoids NaN via the classic Williams "Efficient and Robust Ray-Box"
// trick.
inline bool slabTest(const float bboxMin[3], const float bboxMax[3],
                     const float origin[3], const float invDir[3],
                     float tMax, float& tEntry)
{
    float tx1 = (bboxMin[0] - origin[0]) * invDir[0];
    float tx2 = (bboxMax[0] - origin[0]) * invDir[0];
    float tmin = std::min(tx1, tx2);
    float tmax = std::max(tx1, tx2);

    float ty1 = (bboxMin[1] - origin[1]) * invDir[1];
    float ty2 = (bboxMax[1] - origin[1]) * invDir[1];
    tmin = std::max(tmin, std::min(ty1, ty2));
    tmax = std::min(tmax, std::max(ty1, ty2));

    float tz1 = (bboxMin[2] - origin[2]) * invDir[2];
    float tz2 = (bboxMax[2] - origin[2]) * invDir[2];
    tmin = std::max(tmin, std::min(tz1, tz2));
    tmax = std::min(tmax, std::max(tz1, tz2));

    if (tmax < 0.0f) return false;      // box entirely behind ray
    if (tmin > tmax) return false;      // no overlap
    if (tmin > tMax) return false;      // past the current closest hit
    tEntry = std::max(tmin, 0.0f);
    return true;
}

// Per-triangle scratch used during construction.
struct TriScratch {
    float bmin[3];
    float bmax[3];
    float centroid[3];
};

} // namespace

// --- Construction ------------------------------------------------------------

MeshBVH MeshBVH::build(const MeshData& mesh, int leafSize) {
    MeshBVH bvh;
    if (mesh.indices.empty() || mesh.positions.empty()) return bvh;
    if (leafSize < 1) leafSize = 1;

    const size_t triCount = mesh.triangleCount();
    std::vector<TriScratch> scratch(triCount);
    bvh.triIndices_.resize(triCount);
    for (size_t i = 0; i < triCount; ++i) {
        bvh.triIndices_[i] = (uint32_t)i;
        triBounds(mesh, (uint32_t)i,
                  scratch[i].bmin, scratch[i].bmax, scratch[i].centroid);
    }

    // Reserve a generous upper bound for nodes: worst case ~2N for leafSize=1.
    bvh.nodes_.reserve(triCount * 2 + 1);
    bvh.nodes_.push_back({});  // root

    // Iterative build via an explicit work stack so we don't blow the native
    // stack on huge meshes. Each entry is (nodeIndex, firstTri, triCount).
    struct Work { uint32_t node; uint32_t first; uint32_t count; };
    std::vector<Work> stack;
    stack.reserve(64);
    stack.push_back({0, 0, (uint32_t)triCount});

    const float kInf = std::numeric_limits<float>::infinity();

    while (!stack.empty()) {
        Work w = stack.back();
        stack.pop_back();

        // Compute node bbox + centroid bbox for the assigned triangle range.
        float nodeMin[3] = { kInf, kInf, kInf };
        float nodeMax[3] = { -kInf, -kInf, -kInf };
        float centMin[3] = { kInf, kInf, kInf };
        float centMax[3] = { -kInf, -kInf, -kInf };
        for (uint32_t i = 0; i < w.count; ++i) {
            const TriScratch& t = scratch[bvh.triIndices_[w.first + i]];
            expandBBox(nodeMin, nodeMax, t.bmin, t.bmax);
            expandPoint(centMin, centMax, t.centroid);
        }

        Node& node = bvh.nodes_[w.node];
        for (int a = 0; a < 3; ++a) {
            node.bboxMin[a] = nodeMin[a];
            node.bboxMax[a] = nodeMax[a];
        }

        // Leaf cutoff: few enough tris or degenerate centroid bounds
        // (everything sits on top of itself — further splits can't make
        // progress and would infinite-loop).
        bool degenerate =
            (centMax[0] - centMin[0] < 1e-12f) &&
            (centMax[1] - centMin[1] < 1e-12f) &&
            (centMax[2] - centMin[2] < 1e-12f);
        if (w.count <= (uint32_t)leafSize || degenerate) {
            node.leftFirst = w.first;
            node.triCount  = w.count;
            continue;
        }

        // Split axis = longest centroid extent
        int axis = 0;
        float ext0 = centMax[0] - centMin[0];
        float ext1 = centMax[1] - centMin[1];
        float ext2 = centMax[2] - centMin[2];
        if (ext1 > ext0) { axis = 1; ext0 = ext1; }
        if (ext2 > ext0) { axis = 2; }

        // Median split on centroid along `axis`. Partition triIndices_[first..first+count)
        // in place using nth_element so the left half has the smallest
        // centroid coordinates.
        auto begin = bvh.triIndices_.begin() + w.first;
        auto mid   = begin + w.count / 2;
        auto end   = begin + w.count;
        std::nth_element(begin, mid, end,
            [&](uint32_t a, uint32_t b) {
                return scratch[a].centroid[axis] < scratch[b].centroid[axis];
            });

        uint32_t leftCount  = w.count / 2;
        uint32_t rightCount = w.count - leftCount;
        // Safety: if nth_element produced an empty side (shouldn't happen for
        // count >= 2) fall back to a leaf.
        if (leftCount == 0 || rightCount == 0) {
            node.leftFirst = w.first;
            node.triCount  = w.count;
            continue;
        }

        uint32_t leftChild = (uint32_t)bvh.nodes_.size();
        bvh.nodes_.push_back({});
        bvh.nodes_.push_back({});
        // Re-fetch reference to `node` — the push_back calls may have
        // reallocated the vector, invalidating the earlier reference.
        Node& parent = bvh.nodes_[w.node];
        parent.leftFirst = leftChild;
        parent.triCount  = 0;

        stack.push_back({ leftChild + 1, w.first + leftCount, rightCount });
        stack.push_back({ leftChild,     w.first,             leftCount  });
    }

    return bvh;
}

// --- Queries -----------------------------------------------------------------

BBox MeshBVH::bounds() const {
    BBox b;
    if (nodes_.empty()) return b;
    for (int a = 0; a < 3; ++a) {
        b.min[a] = nodes_[0].bboxMin[a];
        b.max[a] = nodes_[0].bboxMax[a];
    }
    return b;
}

RayHit MeshBVH::raycast(const MeshData& mesh,
                        const float* origin, const float* direction,
                        float maxDistance) const
{
    RayHit result;
    if (nodes_.empty() || !origin || !direction) return result;
    if (mesh.empty() || mesh.indices.empty()) return result;

    float invDir[3];
    for (int a = 0; a < 3; ++a) {
        // Avoid 0 — replace with a very large finite value so slabTest's
        // min/max logic still works. (Strictly infinity would also work on
        // IEEE hardware, but this is the portable version.)
        invDir[a] = (std::fabs(direction[a]) > 1e-30f)
            ? 1.0f / direction[a]
            : std::copysign(1e30f, direction[a] >= 0.0f ? 1.0f : -1.0f);
    }

    float bestT = (maxDistance > 0.0f) ? maxDistance : std::numeric_limits<float>::max();
    float bestU = 0, bestV = 0;
    uint32_t bestTri = 0;
    bool found = false;

    // Iterative traversal with a small fixed stack. BVHs past depth ~64 are
    // pathological; we size for that.
    uint32_t stack[128];
    int sp = 0;
    stack[sp++] = 0;

    while (sp > 0) {
        uint32_t nodeIdx = stack[--sp];
        const Node& n = nodes_[nodeIdx];

        float tEntry;
        if (!slabTest(n.bboxMin, n.bboxMax, origin, invDir, bestT, tEntry))
            continue;

        if (n.triCount > 0) {
            // Leaf: test each triangle.
            for (uint32_t i = 0; i < n.triCount; ++i) {
                uint32_t tri = triIndices_[n.leftFirst + i];
                uint32_t i0 = mesh.indices[tri * 3 + 0];
                uint32_t i1 = mesh.indices[tri * 3 + 1];
                uint32_t i2 = mesh.indices[tri * 3 + 2];
                float t, u, v;
                if (rayTriangle(origin, direction,
                                &mesh.positions[i0 * 3],
                                &mesh.positions[i1 * 3],
                                &mesh.positions[i2 * 3],
                                t, u, v) && t < bestT) {
                    bestT = t;
                    bestU = u;
                    bestV = v;
                    bestTri = tri;
                    found = true;
                }
            }
        } else {
            // Internal: push both children. Push the farther one first so the
            // nearer one is traversed first — gives tighter bestT earlier and
            // prunes more aggressively.
            uint32_t leftIdx  = n.leftFirst;
            uint32_t rightIdx = n.leftFirst + 1;
            const Node& L = nodes_[leftIdx];
            const Node& R = nodes_[rightIdx];
            float tL, tR;
            bool hitL = slabTest(L.bboxMin, L.bboxMax, origin, invDir, bestT, tL);
            bool hitR = slabTest(R.bboxMin, R.bboxMax, origin, invDir, bestT, tR);
            if (hitL && hitR) {
                if (tL < tR) {
                    if (sp + 2 <= (int)(sizeof(stack) / sizeof(stack[0]))) {
                        stack[sp++] = rightIdx;
                        stack[sp++] = leftIdx;
                    }
                } else {
                    if (sp + 2 <= (int)(sizeof(stack) / sizeof(stack[0]))) {
                        stack[sp++] = leftIdx;
                        stack[sp++] = rightIdx;
                    }
                }
            } else if (hitL) {
                if (sp < (int)(sizeof(stack) / sizeof(stack[0])))
                    stack[sp++] = leftIdx;
            } else if (hitR) {
                if (sp < (int)(sizeof(stack) / sizeof(stack[0])))
                    stack[sp++] = rightIdx;
            }
        }
    }

    if (!found) return result;

    // Fill in hit result. Layout mirrors bromesh::raycast so callers can swap
    // between BVH and brute-force without changing downstream code.
    result.hit = true;
    result.distance = bestT;
    result.triangleIndex = bestTri;
    result.baryU = 1.0f - bestU - bestV;
    result.baryV = bestU;
    result.baryW = bestV;
    result.position[0] = origin[0] + direction[0] * bestT;
    result.position[1] = origin[1] + direction[1] * bestT;
    result.position[2] = origin[2] + direction[2] * bestT;

    uint32_t i0 = mesh.indices[bestTri * 3 + 0];
    uint32_t i1 = mesh.indices[bestTri * 3 + 1];
    uint32_t i2 = mesh.indices[bestTri * 3 + 2];
    if (mesh.hasNormals()) {
        for (int c = 0; c < 3; ++c) {
            result.normal[c] = result.baryU * mesh.normals[i0 * 3 + c] +
                               result.baryV * mesh.normals[i1 * 3 + c] +
                               result.baryW * mesh.normals[i2 * 3 + c];
        }
        float len = std::sqrt(result.normal[0] * result.normal[0] +
                              result.normal[1] * result.normal[1] +
                              result.normal[2] * result.normal[2]);
        if (len > 1e-8f) {
            result.normal[0] /= len;
            result.normal[1] /= len;
            result.normal[2] /= len;
        }
    } else {
        const float* p0 = &mesh.positions[i0 * 3];
        const float* p1 = &mesh.positions[i1 * 3];
        const float* p2 = &mesh.positions[i2 * 3];
        float e1[3] = { p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] };
        float e2[3] = { p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2] };
        result.normal[0] = e1[1] * e2[2] - e1[2] * e2[1];
        result.normal[1] = e1[2] * e2[0] - e1[0] * e2[2];
        result.normal[2] = e1[0] * e2[1] - e1[1] * e2[0];
        float len = std::sqrt(result.normal[0] * result.normal[0] +
                              result.normal[1] * result.normal[1] +
                              result.normal[2] * result.normal[2]);
        if (len > 1e-8f) {
            result.normal[0] /= len;
            result.normal[1] /= len;
            result.normal[2] /= len;
        }
    }

    return result;
}

bool MeshBVH::raycastTest(const MeshData& mesh,
                          const float* origin, const float* direction,
                          float maxDistance) const
{
    if (nodes_.empty() || !origin || !direction) return false;
    if (mesh.empty() || mesh.indices.empty()) return false;

    float invDir[3];
    for (int a = 0; a < 3; ++a) {
        invDir[a] = (std::fabs(direction[a]) > 1e-30f)
            ? 1.0f / direction[a]
            : std::copysign(1e30f, direction[a] >= 0.0f ? 1.0f : -1.0f);
    }

    float maxT = (maxDistance > 0.0f) ? maxDistance : std::numeric_limits<float>::max();

    uint32_t stack[128];
    int sp = 0;
    stack[sp++] = 0;

    while (sp > 0) {
        uint32_t nodeIdx = stack[--sp];
        const Node& n = nodes_[nodeIdx];

        float tEntry;
        if (!slabTest(n.bboxMin, n.bboxMax, origin, invDir, maxT, tEntry))
            continue;

        if (n.triCount > 0) {
            for (uint32_t i = 0; i < n.triCount; ++i) {
                uint32_t tri = triIndices_[n.leftFirst + i];
                uint32_t i0 = mesh.indices[tri * 3 + 0];
                uint32_t i1 = mesh.indices[tri * 3 + 1];
                uint32_t i2 = mesh.indices[tri * 3 + 2];
                float t, u, v;
                if (rayTriangle(origin, direction,
                                &mesh.positions[i0 * 3],
                                &mesh.positions[i1 * 3],
                                &mesh.positions[i2 * 3],
                                t, u, v) && t < maxT) {
                    return true;
                }
            }
        } else {
            uint32_t leftIdx  = n.leftFirst;
            uint32_t rightIdx = n.leftFirst + 1;
            if (sp + 2 <= (int)(sizeof(stack) / sizeof(stack[0]))) {
                stack[sp++] = rightIdx;
                stack[sp++] = leftIdx;
            }
        }
    }

    return false;
}

} // namespace bromesh
