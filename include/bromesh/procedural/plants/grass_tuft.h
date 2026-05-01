#pragma once

#include "bromesh/procedural/plants/plant_result.h"

#include <cstdint>

namespace bromesh {

struct GrassTuftParams {
    uint64_t seed = 0;
    int bladeCount = 12;
    float height = 0.4f;
    float baseRadius = 0.06f;
    /// Width of each blade at its base.
    float bladeWidth = 0.015f;
    /// Bend angle at the tip, in radians.
    float bend = 0.6f;
};

/// A clump of curved grass blades. Each blade is a tapered sweep of a thin
/// diamond profile along an arc. No leaves output (blades are the geometry).
PlantResult buildGrassTuft(const GrassTuftParams& params);

} // namespace bromesh
