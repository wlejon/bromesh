#pragma once

#include "bromesh/mesh_data.h"
#include "bromesh/procedural/vec_math.h"

#include <vector>

namespace bromesh {

/// Single placement of a leaf billboard or mesh, intended for instanced
/// rendering. The leaf atlas / variant index is the caller's concern.
struct LeafInstance {
    Vec3 position{};
    Quat orientation{};
    float scale = 1.0f;
    int variantIndex = 0;
};

/// Aggregate output of a plant builder: a single trunk+branches mesh and a
/// list of leaf instances, plus the plant's overall axis-aligned bounding
/// box for culling and placement.
struct PlantResult {
    MeshData branchMesh;
    std::vector<LeafInstance> leaves;
    Vec3 aabbMin{};
    Vec3 aabbMax{};
};

} // namespace bromesh
