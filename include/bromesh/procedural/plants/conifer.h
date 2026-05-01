#pragma once

#include "bromesh/procedural/plants/plant_result.h"

#include <cstdint>

namespace bromesh {

struct ConiferParams {
    uint64_t seed = 0;
    float height = 8.0f;
    float baseRadius = 0.22f;
    /// Number of horizontal whorls (rings of branches) along the trunk.
    int whorlCount = 8;
    /// Branches per whorl.
    int branchesPerWhorl = 6;
    float age01 = 1.0f;
};

/// Build a whorled-branching conifer (single straight leader, periodic
/// rings of side branches that taper toward the top).
PlantResult buildConifer(const ConiferParams& params);

} // namespace bromesh
