#pragma once

#include "bromesh/mesh_data.h"
#include <cstdint>

namespace bromesh {

/// Naive surface nets isosurface extraction.
/// Produces smoother results than marching cubes with shared vertices.
///
/// Sign convention: standard SDF — `field < isoLevel` is inside,
/// `field >= isoLevel` is outside. Normals point outward.
///
/// field: scalar values on a regular grid, row-major.
/// activeBlockMask: optional 8x8x8 cell macro-block activity mask (1 = active, 0 = skip).
MeshData surfaceNets(const float* field, int gridX, int gridY, int gridZ,
                     float isoLevel = 0.0f, float cellSize = 1.0f,
                     const uint8_t* activeBlockMask = nullptr);

} // namespace bromesh
