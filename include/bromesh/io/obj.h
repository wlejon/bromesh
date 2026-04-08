#pragma once

#include "bromesh/mesh_data.h"
#include <string>

namespace bromesh {

/// Load a Wavefront OBJ file. Returns empty mesh on failure.
MeshData loadOBJ(const std::string& path);

/// Save a mesh as Wavefront OBJ.
bool saveOBJ(const MeshData& mesh, const std::string& path);

} // namespace bromesh
