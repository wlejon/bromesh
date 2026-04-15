#pragma once

#include "bromesh/mesh_data.h"

#include <vector>

namespace bromesh {

/// Options for bounded biharmonic weights (Jacobson et al. 2011). Opt-in
/// high-quality path. Requires a clean (manifold or near-manifold) mesh
/// and pays a quadratic-programming cost (minutes on ~50k verts × dozens
/// of bones), so this is not the default — see `autoSelectWeightingMethod`.
///
/// Requires the OSQP submodule (BROMESH_HAS_OSQP). Returns an empty
/// SkinData (with matching boneCount / inverseBindMatrices) if OSQP is
/// not available or the solver fails.
struct BBWOptions {
    /// Max non-zero bones per vertex.
    int maxInfluences = 4;
    /// Weights below this are dropped before final renormalization.
    float minWeight = 1e-3f;
    /// Anchors per bone. The K nearest mesh vertices to the bone's line
    /// segment are constrained to w = 1 for that bone (and w = 0 for the
    /// others). More anchors = smoother weights, at the cost of a less
    /// local stencil.
    int anchorsPerBone = 3;
    /// OSQP tolerances.
    double eps = 1e-4;
    int maxIter = 5000;
};

/// Compute bounded biharmonic skin weights. Each bone yields one quadratic
/// program: minimize the biharmonic energy subject to box constraints
/// (0 ≤ w ≤ 1) and Dirichlet anchors (w = 1 on this bone's anchors, w = 0
/// on other bones' anchors). Weights are stacked, top-K pruned, and
/// renormalized per vertex.
///
/// Returns empty SkinData (with inverseBindMatrices populated) if OSQP is
/// unavailable or any solve fails catastrophically. Callers can fall back
/// to a different weighting method.
SkinData bbwWeights(const MeshData& mesh,
                    const Skeleton& skeleton,
                    const BBWOptions& opts = {});

} // namespace bromesh
