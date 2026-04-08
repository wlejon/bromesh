#pragma once

#include "bromesh/mesh_data.h"

#include <cstdint>

namespace bromesh {

/// Result of surface sampling: a point cloud with positions and normals.
/// Uses MeshData for convenience (indices will be empty).
/// Colors are not populated. UVs are populated if the source mesh has UVs.

/// Sample random points uniformly distributed on the mesh surface.
/// Returns a MeshData with positions (and normals/UVs if the source has them).
/// Indices will be empty (point cloud).
/// seed: random seed for reproducibility (0 = non-deterministic).
MeshData sampleSurface(const MeshData& mesh, size_t numSamples,
                        uint32_t seed = 0);

/// Compute the total surface area of the mesh.
float computeSurfaceArea(const MeshData& mesh);

/// Compute per-triangle areas. Returns one float per triangle.
std::vector<float> computeTriangleAreas(const MeshData& mesh);

} // namespace bromesh
