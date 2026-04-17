#pragma once

#include "bromesh/mesh_data.h"

#include <vector>

namespace bromesh {

// ---------------------------------------------------------------------------
// Polygon triangulation
//
// Turn a 2D polygon outline (optionally with holes) into a triangulated
// mesh. Uses manifold's Triangulate() under the hood, so polygons with holes
// and multiple nested contours are supported.
//
// Empty MeshData is returned if:
//   - the input is degenerate / self-intersecting
//   - manifold isn't linked (BROMESH_HAS_MANIFOLD==0)
//   - the outer contour has fewer than 3 vertices
// ---------------------------------------------------------------------------

/// Triangulate a 2D polygon (with optional holes) into a mesh in the XY
/// plane at z=`z`.
///
/// - `outer`: outer contour as flat x,y,x,y,...  Must be simple (no self-
///   intersections). Wind CCW to emit a face whose normal points toward +Z.
/// - `holes`: optional inner contours, each flat x,y,x,y,... Wind CW
///   relative to the outer contour.
/// - `z`: Z coordinate for all emitted vertices (default 0).
///
/// The returned mesh has per-vertex normals set to [0,0,1] (or [0,0,-1] if
/// the outer contour is wound CW). No UVs or colors are written.
MeshData triangulatePolygon2D(
    const std::vector<float>& outer,
    const std::vector<std::vector<float>>& holes = {},
    float z = 0.0f);

/// Triangulate a planar 3D polygon. Points are projected onto a plane
/// defined by `normal` (must be unit-length), triangulated in 2D, then
/// indexed back to the original 3D vertices.
///
/// - `outer`: outer contour as flat x,y,z,x,y,z,... of 3D points. All points
///   must lie on (or near) a common plane with the given normal.
/// - `holes`: optional inner contours as flat x,y,z,... arrays.
/// - `normal`: 3-component unit plane normal. Triangles are wound so their
///   face normal points along this direction when the outer contour is wound
///   CCW as seen from +normal.
///
/// The returned mesh reuses the input 3D positions verbatim (no projection
/// round-trip) and fills normals with `normal`. Caller owns input/output
/// buffers.
MeshData triangulatePolygon3D(
    const std::vector<float>& outer,
    const std::vector<std::vector<float>>& holes,
    const float normal[3]);

} // namespace bromesh
