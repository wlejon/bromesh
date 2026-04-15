#pragma once

#include "bromesh/mesh_data.h"

#include <vector>

namespace bromesh {

/// Options for surface bone-heat weighting. This is Blender's default auto-
/// weighting algorithm (Baran & Popović 2007, heat-equilibrium formulation).
/// Requires a clean (manifold or near-manifold) triangle mesh; the caller is
/// responsible for deciding whether the mesh is suitable — see
/// `autoSelectWeightingMethod` in weighting.h.
struct BoneHeatOptions {
    /// Max non-zero bones per vertex (glTF standard is 4).
    int maxInfluences = 4;
    /// Weights below this are dropped before final renormalization.
    float minWeight = 1e-3f;
    /// Heat-source multiplier. Higher = tighter influence regions. Baran
    /// & Popović suggest ~1.0; Blender's value is in the same ballpark.
    /// Scales relative to a reference inverse-square distance so the
    /// algorithm is robust to mesh scale.
    float heatStrength = 1.0f;
    /// CG solver tolerance for the per-bone linear solves.
    double solverTol = 1e-7;
    /// CG solver iteration cap.
    int solverMaxIter = 2000;
};

/// Compute skin weights by solving a heat-equilibrium PDE on the mesh
/// surface with each bone acting as a Dirichlet-like source. One sparse
/// SPD solve per bone (CG with Jacobi preconditioning); weights are
/// stacked, top-K pruned, and renormalized per vertex.
///
/// Output SkinData carries the skeleton's inverseBindMatrices verbatim so
/// it is directly usable with applySkinning.
SkinData boneHeatWeights(const MeshData& mesh,
                         const Skeleton& skeleton,
                         const BoneHeatOptions& opts = {});

} // namespace bromesh
