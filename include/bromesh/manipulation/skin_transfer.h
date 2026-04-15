#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Transfer skin weights from a source mesh's SkinData to a target mesh by
/// closest-point projection.
///
/// For each vertex v in `target`, finds the closest point on `source`,
/// interpolates the 3 source-triangle vertex weights by barycentric
/// coordinates, keeps the top-4 bones, and normalizes to sum to 1.
///
/// The returned SkinData uses the same bone index space as `sourceSkin` (so
/// the same Skeleton applies). This is the core enabler for swappable armor:
/// any mesh of similar shape to the base body gets auto-skinned to the same
/// skeleton.
///
/// `maxDistance` of 0 means no distance limit. Vertices farther than
/// `maxDistance` from the source get all-zero weights (and should be
/// filtered or fall back to a default bone by the caller).
SkinData transferSkinWeights(const MeshData& target,
                              const MeshData& source,
                              const SkinData& sourceSkin,
                              float maxDistance = 0.0f);

} // namespace bromesh
