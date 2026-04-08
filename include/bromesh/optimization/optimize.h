#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Optimize vertex/index order for GPU vertex cache.
/// Reorders indices for better post-transform cache utilization.
void optimizeVertexCache(MeshData& mesh);

/// Optimize vertex order for vertex fetch.
/// Reorders vertices to match index access pattern.
void optimizeVertexFetch(MeshData& mesh);

/// Optimize index order to reduce overdraw.
/// Requires a view position (camera location).
void optimizeOverdraw(MeshData& mesh, float threshold = 1.05f);

} // namespace bromesh
