#pragma once

#include "bromesh/mesh_data.h"
#include <vector>

namespace bromesh {

/// Split a mesh into connected components.
/// Returns one MeshData per connected component, each with its own indices
/// starting from 0.
std::vector<MeshData> splitConnectedComponents(const MeshData& mesh);

} // namespace bromesh
