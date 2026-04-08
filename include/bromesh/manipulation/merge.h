#pragma once

#include "bromesh/mesh_data.h"
#include <vector>

namespace bromesh {

/// Merge multiple meshes into a single MeshData.
/// Concatenates all vertex attributes and remaps indices.
/// Only includes attribute streams (normals, UVs, colors) if ALL input meshes have them.
MeshData mergeMeshes(const MeshData* meshes, size_t count);

/// Convenience overload for vector of meshes.
MeshData mergeMeshes(const std::vector<MeshData>& meshes);

} // namespace bromesh
