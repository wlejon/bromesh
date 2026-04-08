#pragma once

#include "bromesh/mesh_data.h"
#include <vector>

namespace bromesh {

/// Parameters for convex decomposition.
struct ConvexDecompParams {
    int maxHulls = 16;           // Maximum number of convex hulls
    int maxVerticesPerHull = 64; // Maximum vertices per hull
    float resolution = 100000;   // Voxel resolution for decomposition
    float minVolumePerHull = 0.001f;
};

/// Decompose a mesh into a set of convex hulls.
/// Uses approximate convex decomposition (V-HACD algorithm).
/// Each returned MeshData is a convex hull suitable for physics collision.
std::vector<MeshData> convexDecomposition(const MeshData& mesh,
                                          const ConvexDecompParams& params = {});

/// Compute the convex hull of a mesh (single hull).
MeshData convexHull(const MeshData& mesh);

} // namespace bromesh
