#pragma once

#include "bromesh/mesh_data.h"
#include <string>

namespace bromesh {

/// Save a mesh as binary STL.
bool saveSTL(const MeshData& mesh, const std::string& path);

/// Load a binary STL file. Returns empty mesh on failure.
MeshData loadSTL(const std::string& path);

} // namespace bromesh
