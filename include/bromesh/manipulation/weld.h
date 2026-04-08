#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Weld (merge) vertices that are within epsilon distance of each other.
/// Deduplicates the vertex buffer and remaps indices.
MeshData weldVertices(const MeshData& mesh, float epsilon = 1e-5f);

} // namespace bromesh
