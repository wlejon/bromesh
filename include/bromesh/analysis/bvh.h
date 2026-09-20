#pragma once

#include "bromesh/mesh_data.h"
#include "bromesh/analysis/raycast.h"

#include <cstdint>
#include <type_traits>
#include <vector>

namespace bromesh {

/// Flat, array-based bounding volume hierarchy over a mesh's triangles.
///
/// Accelerates ray-mesh queries (raycast, raycastTest) for meshes with many
/// triangles. Construction is top-down median-split along the longest axis of
/// the centroid bounds; traversal is an iterative stack using slab AABB tests.
///
/// The BVH indexes triangles by their position in MeshData::indices, so it
/// remains valid only as long as the source mesh's index/position arrays are
/// unchanged. Rebuild whenever the mesh geometry changes.
class MeshBVH {
public:
    MeshBVH() = default;

    /// Build a BVH for `mesh`. `leafSize` is the maximum triangle count per
    /// leaf node (tune for speed/memory; 8 is a reasonable default).
    static MeshBVH build(const MeshData& mesh, int leafSize = 8);

    bool empty() const { return nodes_.empty(); }
    size_t nodeCount() const { return nodes_.size(); }
    size_t triangleCount() const { return triIndices_.size(); }

    /// Overall bounding box of the mesh (the root node's AABB). Zero-sized if
    /// the BVH is empty.
    bromath::AABB3 bounds() const;

    /// Cast a ray against the mesh using this BVH. Semantics match
    /// bromesh::raycast: direction need not be normalized; `maxDistance` of 0
    /// means unlimited. `mesh` must be the same MeshData the BVH was built
    /// against.
    RayHit raycast(const MeshData& mesh,
                   const float* origin, const float* direction,
                   float maxDistance = 0.0f) const;

    /// Fast early-out hit test. Returns true on the first triangle hit.
    bool raycastTest(const MeshData& mesh,
                     const float* origin, const float* direction,
                     float maxDistance = 0.0f) const;

    /// Find the closest point on the mesh surface to `point`.
    RayHit closestPoint(const MeshData& mesh, const float* point) const;

    /// Query all triangle indices whose leaf bounding boxes overlap [boxMin, boxMax].
    /// Invokes `fn(uint32_t triIndex)`. Return false from `fn` to abort traversal early.
    /// Returns false if aborted early, true if traversal completed.
    template <typename Fn>
    bool queryBox(const float boxMin[3], const float boxMax[3], Fn&& fn) const {
        if (nodes_.empty()) return true;
        uint32_t stack[128];
        int sp = 0;
        stack[sp++] = 0;
        while (sp > 0) {
            uint32_t nodeIdx = stack[--sp];
            const Node& n = nodes_[nodeIdx];
            if (boxMin[0] > n.bboxMax[0] || boxMax[0] < n.bboxMin[0] ||
                boxMin[1] > n.bboxMax[1] || boxMax[1] < n.bboxMin[1] ||
                boxMin[2] > n.bboxMax[2] || boxMax[2] < n.bboxMin[2]) {
                continue;
            }
            if (n.triCount > 0) {
                for (uint32_t i = 0; i < n.triCount; ++i) {
                    if constexpr (std::is_same_v<std::invoke_result_t<Fn, uint32_t>, void>) {
                        fn(triIndices_[n.leftFirst + i]);
                    } else {
                        if (!fn(triIndices_[n.leftFirst + i])) {
                            return false;
                        }
                    }
                }
            } else {
                if (sp + 2 <= 128) {
                    stack[sp++] = n.leftFirst + 1;
                    stack[sp++] = n.leftFirst;
                }
            }
        }
        return true;
    }

private:
    /// Flat BVH node. Internal nodes have triCount==0 and leftFirst pointing
    /// at the left child (right child is always leftFirst+1 because children
    /// are allocated in contiguous pairs). Leaves have triCount>0 and
    /// leftFirst is the first index into triIndices_.
    struct Node {
        float bboxMin[3];
        float bboxMax[3];
        uint32_t leftFirst = 0;
        uint32_t triCount  = 0;  // 0 = internal node
    };

    std::vector<Node> nodes_;
    std::vector<uint32_t> triIndices_;  // permutation of [0..mesh.triangleCount())
};

} // namespace bromesh
