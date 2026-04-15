#pragma once

#include "bromesh/mesh_data.h"
#include "bromesh/rigging/rig_spec.h"

#include <string>
#include <vector>

namespace bromesh {

/// A kinematic chain from a "hip" (fork point in the skeleton) down to an
/// IK-chain tip bone (`BoneDecl::ikChain == true`). The hip is the lowest
/// ancestor that has more than one child — typically the pelvis, chest,
/// thorax, or body bone.
struct LegChain {
    std::vector<int> bones;   // hip -> ... -> tip, indices into Skeleton::bones
    int tipBone = -1;         // == bones.back()
    bool grounded = false;    // tip sits in the lower half of the skeleton AABB
};

/// Gait definition. `phases[i]` is the phase offset in [0, 1) for the i-th
/// grounded leg, in the order returned by identifyLegChains (grounded-only,
/// sorted front-to-back then left-to-right). `dutyFactor` is the fraction of
/// the cycle the foot is planted on the ground (0 < dutyFactor <= 1).
struct GaitPattern {
    std::string name;
    std::vector<float> phases;
    float dutyFactor = 0.6f;
};

/// Parameters for a single synthesized cycle.
struct LocomotionParams {
    GaitPattern gait;              // if phases is empty, defaultGait(n) is used
    float strideLength      = 0.30f;
    float cycleDuration     = 1.0f;
    float footLiftHeight    = 0.08f;
    int   keyframesPerCycle = 24;
    float bodyBobAmplitude  = 0.02f; // 0 disables the subtle root-bone bob
    /// Peak swing angle (radians) for each non-grounded IK chain (arms). The
    /// chain's root bone rotates about the world side axis with this
    /// amplitude, phased opposite to the same-side leg so the arms
    /// counter-swing. 0 disables.
    float armSwingAmplitude = 0.35f; // ~20°
    float forwardAxis[3]    = {0.0f, 0.0f, 1.0f};
    float upAxis[3]         = {0.0f, 1.0f, 0.0f};
};

/// Walk the skeleton to enumerate every IK chain declared by the spec. Each
/// chain runs from the nearest fork ancestor (hip/shoulder) down to the tip.
/// Legs are classified "grounded" when the tip sits in the lower half of the
/// skeleton's bind-pose AABB along the up axis. The returned vector contains
/// *all* chains; grounded-only callers should filter.
std::vector<LegChain> identifyLegChains(const Skeleton& skeleton,
                                        const RigSpec& spec,
                                        const float upAxis[3] = nullptr);

/// Built-in gaits for the common leg counts:
///   2 -> biped walk      (alternating)
///   4 -> quadruped trot  (diagonal pairs)
///   6 -> hexapod tripod  (alternating triangles)
///   8 -> octopod ripple  (front-to-back wave, L/R anti-phase)
/// Any other leg count returns an empty GaitPattern — callers must supply a
/// custom gait.
GaitPattern defaultGait(int groundedLegCount);

/// Generate a looping walk/run cycle for the given skeleton. The returned
/// Animation is the same type produced by loadGLTF — ready to drive through
/// evaluateAnimation, blendPoses, saveGLTF, etc.
///
/// If `params.gait.phases` is empty, defaultGait(groundedLegCount) is
/// substituted. If no gait can be determined (unusual leg counts, no
/// grounded legs), an empty Animation is returned.
Animation generateLocomotionCycle(const Skeleton& skeleton,
                                  const RigSpec& spec,
                                  const LocomotionParams& params);

} // namespace bromesh
