#pragma once

#include "bromesh/procedural/plants/plant_result.h"

#include <cstdint>

namespace bromesh {

struct SucculentParams {
    uint64_t seed = 0;
    /// Number of leaves in the rosette.
    int leafCount = 24;
    /// Maximum leaf length.
    float leafLength = 0.35f;
    /// Maximum half-width of a leaf at its widest point.
    float leafWidth = 0.07f;
    /// Maximum thickness (Y axis of profile).
    float leafThickness = 0.04f;
    /// Tilt of leaves away from horizontal, radians.
    float tilt = 0.6f;
};

/// Rosette-form succulent: thick tapered leaves arranged in a phyllotactic
/// spiral around a central origin. Pure swept geometry, no leaf instances.
PlantResult buildSucculent(const SucculentParams& params);

} // namespace bromesh
