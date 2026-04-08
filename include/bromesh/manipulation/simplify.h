#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Simplify a mesh using quadric error metrics.
/// targetRatio: fraction of triangles to keep (0.0 to 1.0).
/// targetError: maximum allowed error (0 = no limit).
/// Returns a new simplified mesh.
MeshData simplify(const MeshData& mesh, float targetRatio,
                  float targetError = 0.01f);

/// Simplify a mesh using quadric error metrics with attribute preservation.
/// Incorporates UV and normal attributes into the error metric so that
/// UV seams and hard edges are preserved during simplification.
/// uvWeight/normalWeight control relative importance (0 = ignore).
MeshData simplifyWithAttributes(const MeshData& mesh, float targetRatio,
                                float targetError = 0.01f,
                                float uvWeight = 1.0f,
                                float normalWeight = 0.5f);

/// Generate an LOD chain: multiple simplification levels.
/// ratios: array of target ratios (e.g. {0.5, 0.25, 0.125}).
/// Returns one MeshData per ratio.
std::vector<MeshData> generateLODChain(const MeshData& mesh,
                                       const float* ratios, int count);

} // namespace bromesh
