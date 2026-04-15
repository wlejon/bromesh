#pragma once

#include "bromesh/mesh_data.h"

#include <cstdint>
#include <vector>

namespace bromesh {

/// A pose is a flat array of local transforms, one per bone, in the same
/// index space as Skeleton::bones:
///   translation: 3 floats (offset 0)
///   rotation:    4 floats (quaternion xyzw, offset 3)
///   scale:       3 floats (offset 7)
/// Stride = 10 floats per bone.
struct Pose {
    std::vector<float> data;
    size_t boneCount() const { return data.size() / 10; }
};

/// Initialize a pose from the skeleton's bind (rest) transforms.
Pose bindPose(const Skeleton& skeleton);

/// Evaluate an animation at time `t` (seconds). Starts from `skeleton` bind
/// pose, then overrides only the bones that have channels in the animation.
/// Time clamps to [0, anim.duration]; use `loop=true` to wrap instead.
Pose evaluateAnimation(const Skeleton& skeleton,
                       const Animation& anim,
                       float t,
                       bool loop = true);

/// In-place variant. `pose` must already be sized to skeleton.bones.size().
void evaluateAnimationInto(const Skeleton& skeleton,
                           const Animation& anim,
                           float t,
                           bool loop,
                           Pose& pose);

/// Blend `b` into `a` by weight in [0,1] (0 = all a, 1 = all b).
/// If `boneMask` is non-null it must be length == boneCount; a 1 enables
/// blending on that bone, a 0 keeps `a`'s value. Useful for upper-body vs
/// lower-body layering.
/// Rotations blend as slerp, translations and scales as lerp.
void blendPoses(Pose& a, const Pose& b, float weight,
                const uint8_t* boneMask = nullptr);

/// Compose a pose into world-space bone matrices by walking the skeleton
/// hierarchy. Output stride: 16 floats per bone, column-major.
/// `outWorld` is resized to skeleton.bones.size() * 16.
void computeWorldMatrices(const Skeleton& skeleton,
                          const Pose& pose,
                          std::vector<float>& outWorld);

/// Compose final skinning (joint) matrices: world(bone) * inverseBind(bone).
/// These are what applySkinning consumes. Output is skeleton.bones.size() *
/// 16 floats, column-major, ready for the pose-matrix argument of
/// applySkinning.
void computeSkinningMatrices(const Skeleton& skeleton,
                             const Pose& pose,
                             std::vector<float>& outSkinning);

/// Compute the world-space transform of a named socket given a pose.
/// Returns false if not found. `outMatrix` is 16 floats column-major.
bool socketWorldMatrix(const Skeleton& skeleton,
                       const Pose& pose,
                       const std::string& socketName,
                       float* outMatrix);

} // namespace bromesh
