#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Tunables for weight post-processing.
struct WeightPostProcessOptions {
    /// Laplacian smoothing passes. 0 disables smoothing (outlier rejection
    /// and re-normalization still run).
    int iterations = 4;
    /// Per-iteration blend factor toward the neighbor average. 0 = no change,
    /// 1 = replace each weight with the neighborhood mean. 0.5 is a good
    /// default — removes voxel-grid stairstepping without over-blurring joint
    /// boundaries.
    float alpha = 0.5f;
    /// Weights below this are dropped before final normalization.
    float minWeight = 1e-3f;
};

/// Post-process voxel-bind output to remove stairstepping and spurious
/// per-vertex outliers. Expands the top-4 representation to dense
/// per-vertex × per-bone, runs N passes of Laplacian smoothing over mesh
/// adjacency, then re-selects the top-4 influences and re-normalizes.
///
/// Stable, deterministic, operates in place. Safe to call with iterations=0
/// to get pure re-normalization + min-weight pruning.
void postProcessWeights(const MeshData& mesh,
                        SkinData& skin,
                        const WeightPostProcessOptions& opts = {});

} // namespace bromesh
