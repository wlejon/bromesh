#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Remove degenerate triangles (zero-area: two or more identical vertex indices,
/// or collinear vertices within epsilon).
/// Returns a new mesh with degenerate triangles removed.
MeshData removeDegenerateTriangles(const MeshData& mesh, float areaEpsilon = 1e-8f);

/// Remove duplicate triangles (same three vertex positions regardless of winding).
/// Returns a new mesh with duplicates removed.
MeshData removeDuplicateTriangles(const MeshData& mesh);

/// Fill simple boundary loops (holes) by connecting boundary edges with a fan.
/// Only fills holes with up to maxEdges boundary edges.
/// Returns a new mesh with holes filled. Normals are computed for new triangles.
MeshData fillHoles(const MeshData& mesh, int maxEdges = 64);

} // namespace bromesh
