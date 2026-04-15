#pragma once

#include "bromesh/mesh_data.h"

#include <vector>

namespace bromesh {

/// Tunables for geodesic voxel binding.
struct VoxelBindOptions {
    /// Longest-axis voxel count. Grid resolution on the other axes scales
    /// proportionally to the mesh bounding box.
    int maxResolution = 96;
    /// Max non-zero bones per vertex. glTF standard is 4.
    int maxInfluences = 4;
    /// Falloff exponent for 1 / (dist^p + eps) weighting. Higher = tighter
    /// bone influence regions.
    float falloffPower = 4.0f;
    /// Weights below this are dropped before final renormalization.
    float minWeight = 1e-3f;
};

/// Compute skin weights by geodesic voxel binding (Dionne & de Lasa 2013,
/// simplified). Steps: voxelize the mesh into a solid grid, run a BFS from
/// each bone through the solid voxels to get per-voxel distances, convert
/// to weights, sample per vertex, keep top N influences.
///
/// Output's inverseBindMatrices are copied from `skeleton` so the returned
/// SkinData is directly usable with applySkinning.
SkinData voxelBindWeights(const MeshData& mesh,
                          const Skeleton& skeleton,
                          const VoxelBindOptions& opts = {});

} // namespace bromesh
