#pragma once

#include "bromesh/mesh_data.h"

#include <string>
#include <vector>

namespace bromesh {

/// Remap an animation from `srcSkeleton` to `dstSkeleton` by matching bone
/// names. Channels whose bone name has no counterpart in the destination
/// skeleton are dropped. This is the basic "play one walk cycle on many
/// characters" retarget and assumes both skeletons share the same naming
/// convention (Rigify-like, Mixamo, etc).
///
/// Only index remapping is done — bone rest poses are NOT adjusted. For
/// characters whose rest poses differ significantly, a full retarget
/// pipeline is needed (bone-space delta conversion), which is out of scope
/// here. In practice, Rigify-derived skeletons share rest poses closely
/// enough that name-based remap is a usable first pass.
Animation retargetAnimation(const Animation& anim,
                             const Skeleton& srcSkeleton,
                             const Skeleton& dstSkeleton);

/// Common socket names that this helper will populate from a Rigify-style
/// skeleton: "hand.L", "hand.R", "foot.L", "foot.R", "head", "chest",
/// "pelvis". Looks for the bone whose name exactly matches, with these
/// fallback patterns tried in order:
///   "DEF-<name>" (Rigify deform prefix)
///   "<name>"
///   "mixamorig:<name>" (Mixamo prefix)
///
/// Sockets are appended with identity offset. Returns the number added.
int addRigifySockets(Skeleton& skeleton);

/// Find a bone by name suffix (case-sensitive). Useful when skeletons are
/// exported with a prefix like "DEF-" or "mixamorig:". Returns the first
/// matching bone, or -1 if none.
int findBoneBySuffix(const Skeleton& skeleton, const std::string& suffix);

} // namespace bromesh
