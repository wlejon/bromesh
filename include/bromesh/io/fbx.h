#pragma once

#include "bromesh/mesh_data.h"
#include <string>
#include <vector>

namespace bromesh {

/// Load an FBX file. Returns all meshes found in the scene.
/// Supports positions, normals, UVs, and vertex colors.
/// Returns empty vector on failure.
std::vector<MeshData> loadFBX(const std::string& path);

} // namespace bromesh
