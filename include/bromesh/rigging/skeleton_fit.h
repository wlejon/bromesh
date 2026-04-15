#pragma once

#include "bromesh/mesh_data.h"
#include "bromesh/rigging/landmarks.h"
#include "bromesh/rigging/rig_spec.h"

#include <string>
#include <vector>

namespace bromesh {

/// Fit a Skeleton to `mesh` according to `spec` using `landmarks` to resolve
/// each bone's rest position. `mesh` is currently unused but passed in so
/// future passes (inward projection, surface-following spine, validation)
/// can consume it without an API change.
///
/// Bones in the returned Skeleton are topologically sorted (parents precede
/// children). Each bone's local TRS and inverseBind are set so the bind-pose
/// world transform places +Y along the bone's head→tail direction.
///
/// Landmarks missing from `landmarks` cause the corresponding bones to fall
/// back to a zero-length placeholder; use `missingLandmarks()` before calling
/// if you want a hard failure on missing input.
Skeleton fitSkeleton(const RigSpec& spec,
                     const Landmarks& landmarks,
                     const MeshData& mesh);

} // namespace bromesh
