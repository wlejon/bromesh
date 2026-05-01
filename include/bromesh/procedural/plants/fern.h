#pragma once

#include "bromesh/procedural/plants/plant_result.h"

#include <cstdint>

namespace bromesh {

struct FernParams {
    uint64_t seed = 0;
    /// Number of pinnate leaflet pairs along the frond.
    int leafletPairs = 14;
    /// Frond length.
    float length = 1.5f;
    /// Stem radius at the base.
    float stemRadius = 0.012f;
    /// Maximum leaflet length (mid-frond).
    float leafletLength = 0.35f;
    /// Curvature of the frond (radians of arc total).
    float curvature = 0.8f;
    /// 0..1 frond length fraction. Pinnae count is fixed by leafletPairs;
    /// age trims the rachis from the tip.
    float age01 = 1.0f;
};

/// Single fern frond: an L-system drives a curved central rachis with
/// pairs of pinnate leaflets that taper toward the tip. Demonstrates the
/// LSystem class composed with sweep + leaf instancing.
PlantResult buildFern(const FernParams& params);

} // namespace bromesh
