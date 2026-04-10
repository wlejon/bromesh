#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Classic marching cubes isosurface extraction.
/// field: scalar values on a regular grid, row-major (x varies fastest).
/// gridX/Y/Z: grid dimensions.
/// isoLevel: threshold value for the surface.
/// cellSize: world-space size of each grid cell (default 1.0).
/// closeBoundary: when true (default), virtually pad the field with a
///   sub-iso sentinel layer so the iso-surface is closed at the grid edge.
///   Required for downstream analysis (volume, raycast, BVH) and for
///   terrain chunks where the surface naturally extends past the chunk
///   boundary. Set to false to keep the surface open at the grid edge —
///   useful when stitching chunks together explicitly.
MeshData marchingCubes(const float* field, int gridX, int gridY, int gridZ,
                       float isoLevel = 0.0f, float cellSize = 1.0f,
                       bool closeBoundary = true);

} // namespace bromesh
