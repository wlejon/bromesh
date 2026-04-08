#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Loop subdivision: splits each triangle into 4 using midpoint insertion
/// and smoothing. Works on triangle meshes. Preserves UVs by linear
/// interpolation along edges.
/// iterations: number of subdivision passes (each pass 4x triangle count).
MeshData subdivideLoop(const MeshData& mesh, int iterations = 1);

/// Catmull-Clark subdivision: inserts face and edge points, producing a
/// quad-dominant mesh (output is still triangulated for MeshData compat).
/// iterations: number of subdivision passes.
MeshData subdivideCatmullClark(const MeshData& mesh, int iterations = 1);

/// Simple midpoint subdivision: splits each triangle into 4 without
/// smoothing vertex positions. Pure geometric refinement.
MeshData subdivideMidpoint(const MeshData& mesh, int iterations = 1);

} // namespace bromesh
