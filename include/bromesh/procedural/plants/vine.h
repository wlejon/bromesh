#pragma once

#include "bromesh/procedural/plants/plant_result.h"

#include <cstdint>

namespace bromesh {

struct VineParams {
    uint64_t seed = 0;
    /// Length along the helical path.
    float length = 6.0f;
    /// Vine stem radius.
    float radius = 0.04f;
    /// Helix winding radius.
    float helixRadius = 0.5f;
    /// Number of full turns over the length.
    float turns = 3.0f;
    /// Leaves per unit length.
    float leafDensity = 4.0f;
    /// 0..1 uniform scale of length. The helix structure (turns, helixRadius)
    /// is fixed; age just truncates the path.
    float age01 = 1.0f;
};

/// A long sweeping vine: noisy helical path, periodic leaves placed at
/// regular path arc-length intervals.
PlantResult buildVine(const VineParams& params);

} // namespace bromesh
