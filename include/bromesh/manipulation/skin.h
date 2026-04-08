#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Apply skeletal skinning to a mesh using the given pose matrices.
/// poseMatrices: array of 4x4 column-major matrices, one per bone (boneCount * 16 floats).
/// Transforms positions and normals in-place according to bone weights/indices.
void applySkinning(MeshData& mesh, const SkinData& skin,
                   const float* poseMatrices);

/// Apply a morph target to a mesh at the given weight (0 = no effect, 1 = full).
/// Adds deltaPositions (and deltaNormals if present) scaled by weight.
void applyMorphTarget(MeshData& mesh, const MorphTarget& morph, float weight);

/// Normalize bone weights so they sum to 1.0 for each vertex.
/// Discards near-zero weights and repacks indices.
void normalizeWeights(SkinData& skin);

} // namespace bromesh
