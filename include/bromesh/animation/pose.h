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

/// Weighted N-way pose blend into `out`.
///
/// `poses`/`weights` are parallel arrays of `count` sources. Weights are
/// normalized internally (any positive sum works); if the sum is ~0 the
/// first pose wins. All poses must share one bone count.
///
/// Semantics (the industry-standard N-way accumulate):
///   - translations and scales: weighted sum of the normalized weights.
///   - rotations, count == 2: exact slerp — bit-identical to
///     blendPoses(a, b, w1/(w0+w1)), so two-way blend spaces and the
///     existing crossfade path agree exactly.
///   - rotations, count >= 3: weighted nlerp — each quaternion is
///     hemisphere-aligned against the highest-weight source (negated when
///     its dot with that reference is negative), accumulated by weight,
///     then normalized. Nlerp is chosen over iterative slerp because it is
///     commutative across sources, allocation-free, and matches what game
///     runtimes ship for N-way blending; it coincides with slerp at the
///     endpoints and at 50/50, and stays within a fraction of a degree of
///     it for the < 90 degree bone deltas typical between locomotion clips.
///
/// `out` is resized if needed and initialized from the highest-weight
/// source, so with a `boneMask` (same contract as blendPoses: 1 = blend
/// this bone) the masked-OUT bones keep whatever `out` already held when
/// it was correctly sized, and take the reference pose's value otherwise.
/// `out` must not alias any input pose. No allocation when `out` is
/// already sized.
void blendPosesN(const Pose* const* poses, const float* weights, size_t count,
                 Pose& out, const uint8_t* boneMask = nullptr);

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
