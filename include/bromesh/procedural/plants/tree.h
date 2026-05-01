#pragma once

#include "bromesh/procedural/plants/plant_result.h"

#include <cstdint>

namespace bromesh {

struct TreeParams {
    uint64_t seed = 0;
    /// Maximum height of the canopy from the base of the trunk.
    float height = 6.0f;
    /// Trunk radius at the base.
    float trunkRadius = 0.18f;
    /// Lateral radius of the attractor cloud (canopy reach).
    float canopyRadius = 3.0f;
    /// Number of attractors sampled in the canopy. More -> denser branching.
    int attractorCount = 600;
    /// 0..1 scales height, canopy, and leaf density. Useful for growth LOD.
    float age01 = 1.0f;
};

/// Build a broadleaf tree using space colonization + parallel-transport
/// sweep for branch geometry. Leaves are placed on terminal segments.
PlantResult buildTree(const TreeParams& params);

} // namespace bromesh
