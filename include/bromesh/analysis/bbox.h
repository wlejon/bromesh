#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Compute the axis-aligned bounding box of a mesh.
BBox computeBBox(const MeshData& mesh);

/// Check if a mesh is manifold (watertight, every edge shared by exactly 2 triangles).
bool isManifold(const MeshData& mesh);

/// Compute the volume of a closed (manifold) mesh using the divergence theorem.
/// Returns 0 if the mesh is not closed.
float computeVolume(const MeshData& mesh);

} // namespace bromesh
