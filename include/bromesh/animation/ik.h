#pragma once

#include "bromesh/mesh_data.h"
#include "bromesh/animation/pose.h"

#include <vector>

namespace bromesh {

/// Two-bone analytic IK. Classic "elbow / knee" solver: given a 3-bone chain
/// root -> mid -> end (shoulder -> elbow -> wrist), rotates root and mid so
/// that the end bone's tip reaches `targetWorld` as closely as possible. The
/// end bone's rotation is also aligned toward the target.
///
/// `poleWorld` (optional, may be null) disambiguates the swing plane: the
/// elbow will be pushed toward the pole position. Without it, the mid bone
/// stays in its current plane.
///
/// The solver mutates `pose` in-place. Returns false if the chain indices
/// are invalid or bone lengths are degenerate.
///
/// Bone indices are in Skeleton::bones order. `root` must be an ancestor of
/// `mid`, which must be the parent of `end`. The chain must be rigid (no
/// intermediate bones).
bool solveTwoBoneIK(const Skeleton& skeleton,
                    Pose& pose,
                    int rootBone, int midBone, int endBone,
                    const float targetWorld[3],
                    const float* poleWorld = nullptr);

/// FABRIK (Forward And Backward Reaching IK) for chains of arbitrary length.
/// `chain` is bone indices from root to tip. The tip converges to
/// `targetWorld` via iterative two-pass reaching. Fast, stable, works well
/// for spines, tails, and tentacles.
///
/// `iterations` and `tolerance` are convergence controls; defaults are fine
/// for short chains.
bool solveFABRIK(const Skeleton& skeleton,
                 Pose& pose,
                 const std::vector<int>& chain,
                 const float targetWorld[3],
                 int iterations = 10,
                 float tolerance = 1e-3f);

/// Orient `bone` so that its local +axis points at `targetWorld`. `upAxis`
/// keeps the roll stable. Both axes are in the bone's local space. Useful
/// for head look-at, turret aiming, etc.
bool solveLookAt(const Skeleton& skeleton,
                 Pose& pose,
                 int bone,
                 const float targetWorld[3],
                 const float localForward[3] = nullptr,
                 const float localUp[3] = nullptr);

} // namespace bromesh
