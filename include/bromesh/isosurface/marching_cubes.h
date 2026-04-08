#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Classic marching cubes isosurface extraction.
/// field: scalar values on a regular grid, row-major (x varies fastest).
/// gridX/Y/Z: grid dimensions.
/// isoLevel: threshold value for the surface.
/// cellSize: world-space size of each grid cell (default 1.0).
MeshData marchingCubes(const float* field, int gridX, int gridY, int gridZ,
                       float isoLevel = 0.0f, float cellSize = 1.0f);

} // namespace bromesh
