#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Reorder triangles for spatial locality.
/// Improves performance of subsequent vertex cache optimization.
void spatialSortTriangles(MeshData& mesh);

/// Reorder vertices for spatial locality using a space-filling curve.
/// Improves compression and streaming performance.
void spatialSortVertices(MeshData& mesh);

/// Generate a shadow index buffer optimized for position-only rendering (Z-prepass / shadow maps).
/// Returns an index buffer where all vertices with identical positions share the same index,
/// enabling more efficient rendering when only positions are needed.
std::vector<uint32_t> generateShadowIndexBuffer(const MeshData& mesh);

} // namespace bromesh
