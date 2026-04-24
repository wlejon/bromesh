#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Generate a box mesh centered at origin.
MeshData box(float halfW, float halfH, float halfD);

/// Generate a UV sphere mesh centered at origin.
MeshData sphere(float radius, int segments = 16, int rings = 12);

/// Generate a cylinder mesh along the Y axis, centered at origin.
MeshData cylinder(float radius, float halfHeight, int segments = 16);

/// Generate a capsule mesh along the Y axis, centered at origin.
MeshData capsule(float radius, float halfHeight, int segments = 16, int rings = 8);

/// Generate a flat plane in the XZ plane centered at origin.
MeshData plane(float halfW, float halfD, int subdivX = 1, int subdivZ = 1);

/// Generate a torus mesh centered at origin, lying in the XZ plane.
MeshData torus(float majorRadius, float minorRadius,
               int majorSegments = 24, int minorSegments = 12);

/// Generate a mesh from a heightmap.
/// heights: row-major grid of Y values.
/// gridW/gridH: grid dimensions of the meshed area.
/// cellSize: horizontal spacing between grid points.
/// border: if > 0, `heights` is expected to be padded with `border` extra
///   rows/columns on every side (full array size
///   `(gridW + 2*border) * (gridH + 2*border)`). The outer ring is used only
///   to compute proper central-difference normals at the interior boundary,
///   so chunked terrain can produce matching normals at shared edges. The
///   mesh geometry still covers only the inner `gridW * gridH`.
MeshData heightmapGrid(const float* heights, int gridW, int gridH,
                       float cellSize = 1.0f, int border = 0);

} // namespace bromesh
