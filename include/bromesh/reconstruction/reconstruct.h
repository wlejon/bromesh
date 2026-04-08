#pragma once

#include "bromesh/mesh_data.h"
#include <vector>

namespace bromesh {

/// Parameters for point cloud reconstruction.
struct ReconstructParams {
    int gridResolution = 64;    // Voxel grid resolution along longest axis
    float supportRadius = 0.0f; // Influence radius per point (0 = auto from density)
    float isoLevel = 0.5f;      // Threshold for surface extraction
};

/// Reconstruct a triangle mesh from an oriented point cloud using implicit
/// surface reconstruction. Points must have normals.
/// positions: xyz, stride 3 (pointCount points)
/// normals: xyz, stride 3 (pointCount normals, must be unit length)
/// Returns a triangle mesh via marching cubes on the reconstructed field.
MeshData reconstructFromPointCloud(
    const float* positions, const float* normals, size_t pointCount,
    const ReconstructParams& params = {});

/// Convenience: reconstruct from a MeshData used as point cloud
/// (only positions and normals are used, indices are ignored).
MeshData reconstructFromPointCloud(
    const MeshData& pointCloud,
    const ReconstructParams& params = {});

} // namespace bromesh
