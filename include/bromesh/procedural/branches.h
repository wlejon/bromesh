#pragma once

#include "bromesh/mesh_data.h"
#include "bromesh/procedural/space_colonization.h"

#include <vector>

namespace bromesh {

/// Mesh a branch tree (any segment list with `parent` indices) as one merged
/// mesh. Each maximal chain of single-child segments is meshed as one
/// multi-ring sweep so consecutive rings share parallel-transport orientation
/// — no per-segment ring discontinuities. Forks blend via the parent's
/// continuing radius rather than capping with a flat disc.
///
/// `sides` is the number of cross-sectional segments per ring (e.g. 6 or 8).
MeshData meshBranches(const std::vector<BranchSegment>& segments, int sides);

} // namespace bromesh
