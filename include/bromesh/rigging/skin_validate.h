#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Per-vertex quality metrics for a SkinData. Produced by validateSkin.
///
/// Definitions:
///  - orphan: vertex with no non-zero weight.
///  - badSum: |sum(weights) - 1| > sumTolerance.
///  - nan:    any weight is NaN or negative.
///  - maxSumDeviation: the largest |sum - 1| seen across all vertices.
///  - maxInfluencesObserved: max number of non-zero weights on any vertex.
struct SkinValidation {
    size_t vertexCount = 0;
    size_t orphanCount = 0;
    size_t badSumCount = 0;
    size_t nanCount    = 0;
    float  maxSumDeviation = 0.0f;
    int    maxInfluencesObserved = 0;

    /// True iff no orphan, bad-sum, or NaN vertices were observed.
    bool clean() const {
        return orphanCount == 0 && badSumCount == 0 && nanCount == 0;
    }
};

/// Walk a SkinData and report quality metrics. `influences` is the stride
/// of boneWeights/boneIndices (glTF standard is 4). `sumTolerance` controls
/// what counts as a bad sum — 1e-3 is the bar the auto-rig tests use.
SkinValidation validateSkin(const MeshData& mesh,
                            const SkinData& skin,
                            int influences = 4,
                            float sumTolerance = 1e-3f);

} // namespace bromesh
