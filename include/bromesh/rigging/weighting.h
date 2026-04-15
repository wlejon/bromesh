#pragma once

#include "bromesh/mesh_data.h"
#include "bromesh/rigging/bbw.h"
#include "bromesh/rigging/bone_heat.h"
#include "bromesh/rigging/voxel_bind.h"

namespace bromesh {

/// Which weighting algorithm to use. `Auto` picks at runtime based on mesh
/// manifoldness: manifold → BoneHeat (fast, smooth), non-manifold → VoxelBind
/// (topology-robust). BBW is opt-in only; it's higher quality on clean meshes
/// but costs a QP solve per bone.
enum class WeightingMethod {
    Auto,
    VoxelBind,
    BoneHeat,
    BBW,
};

/// Unified options bag for `autoRig`. Each algorithm-specific block is
/// consulted only when its method is selected.
struct WeightingOptions {
    WeightingMethod method = WeightingMethod::Auto;
    VoxelBindOptions voxel;
    BoneHeatOptions  boneHeat;
    BBWOptions       bbw;
    /// Laplacian weight-smoothing passes applied to the final skin. Applies
    /// to all methods (including BBW — cheap, removes any residual noise).
    int smoothIterations = 2;
    float smoothAlpha = 0.5f;
    /// Minimum weight used by the top-K pruner (and matches minWeight in
    /// the per-method options block). Keep consistent with method.minWeight.
    float minWeight = 1e-3f;
};

/// Pick an appropriate concrete method when the user asks for `Auto`. Checks
/// mesh manifoldness. Returns BoneHeat for clean manifold meshes, VoxelBind
/// otherwise. Never returns Auto / BBW.
WeightingMethod autoSelectWeightingMethod(const MeshData& mesh);

/// Stable string name for a WeightingMethod. Round-trips through
/// parseWeightingMethod. Used at language/data boundaries (JSON, JS bindings).
const char* weightingMethodName(WeightingMethod m);

/// Parse a stable string back to a WeightingMethod. Accepts the exact
/// lowercase names emitted by weightingMethodName. Unknown strings fall
/// back to WeightingMethod::Auto.
WeightingMethod parseWeightingMethod(const char* name);

} // namespace bromesh
