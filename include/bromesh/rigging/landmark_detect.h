#pragma once

#include "bromesh/mesh_data.h"
#include "bromesh/rigging/landmarks.h"

namespace bromesh {

/// Options for geometric landmark detection. Defaults assume the mesh is in
/// T-pose with +Y up and +Z forward. Non-T-pose meshes may still work if
/// limbs are roughly extended, but quality degrades.
struct LandmarkDetectOptions {
    float upAxis[3]      = { 0.0f, 1.0f, 0.0f };
    float forwardAxis[3] = { 0.0f, 0.0f, 1.0f };
    /// Fraction of body height used as foot length when projecting toe
    /// landmarks forward from the ankle.
    float footLengthFrac = 0.08f;
};

/// Derive the 18 landmarks declared by builtinHumanoidSpec() from mesh shape
/// alone. Uses AABB extrema and bilateral symmetry rather than a learned
/// model — good for T-pose humanoids, falls over on highly stylised poses.
///
/// The output is consumable by autoRig() with the humanoid spec and no other
/// work. Use this when no user-marked landmarks are available.
Landmarks detectHumanoidLandmarks(const MeshData& mesh,
                                  const LandmarkDetectOptions& opts = {});

} // namespace bromesh
