#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Dual contouring isosurface extraction.
/// Preserves sharp features using QEF (quadric error function) minimization.
/// field: scalar values on a regular grid, row-major.
MeshData dualContour(const float* field, int gridX, int gridY, int gridZ,
                     float isoLevel = 0.0f, float cellSize = 1.0f);

} // namespace bromesh
