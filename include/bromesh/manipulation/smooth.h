#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Laplacian smoothing: moves each vertex toward the centroid of its neighbors.
/// lambda: smoothing factor per iteration (0.0 to 1.0, typically 0.5).
/// iterations: number of smoothing passes.
/// Modifies vertex positions in-place. Preserves topology.
void smoothLaplacian(MeshData& mesh, float lambda = 0.5f, int iterations = 1);

/// Taubin smoothing: alternating positive/negative Laplacian to prevent shrinkage.
/// lambda: positive smoothing factor (typically 0.5).
/// mu: negative smoothing factor (typically -0.53, must be < -lambda).
/// iterations: number of smooth+inflate cycles.
void smoothTaubin(MeshData& mesh, float lambda = 0.5f, float mu = -0.53f,
                  int iterations = 1);

} // namespace bromesh
