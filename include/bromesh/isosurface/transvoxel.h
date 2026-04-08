#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Transvoxel algorithm for seamless LOD transitions between chunks.
/// Generates transition cells along chunk boundaries where neighboring
/// chunks have different resolution levels.
///
/// field: scalar values for this chunk's grid, row-major.
/// gridSize: cubic grid dimension (must be uniform).
/// lod: LOD level of this chunk (0 = highest detail).
/// neighborLods: LOD level of each neighbor [+X, -X, +Y, -Y, +Z, -Z].
///               -1 means no neighbor (chunk boundary = world edge).
MeshData transvoxel(const float* field, int gridSize, int lod,
                    const int neighborLods[6],
                    float isoLevel = 0.0f, float cellSize = 1.0f);

} // namespace bromesh
