#pragma once

#include "bromesh/mesh_data.h"
#include <string>

namespace bromesh {

/// Save a mesh as binary PLY.
/// Writes positions, normals (if present), UVs (if present), and colors (if present).
bool savePLY(const MeshData& mesh, const std::string& path);

/// Load a PLY file (binary little-endian or ASCII).
/// Returns empty mesh on failure.
MeshData loadPLY(const std::string& path);

} // namespace bromesh
