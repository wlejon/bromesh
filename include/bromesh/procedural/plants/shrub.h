#pragma once

#include "bromesh/procedural/plants/plant_result.h"

#include <cstdint>

namespace bromesh {

struct ShrubParams {
    uint64_t seed = 0;
    float height = 1.5f;
    float radius = 1.2f;
    int stemCount = 4;
    int attractorCount = 200;
    float age01 = 1.0f;
};

/// Multi-stem shrub: several seed points clustered at the base, dense
/// space-colonization growth, leaves on terminal segments.
PlantResult buildShrub(const ShrubParams& params);

} // namespace bromesh
