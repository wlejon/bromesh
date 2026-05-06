#pragma once

#include "bromesh/mesh_data.h"
#include "bromesh/procedural/space_colonization.h"

#include <cstdint>
#include <vector>

namespace bromesh {

/// Options for `placeLeavesOnBranches` / `scatterLeaves`.
struct LeafPlacementOptions {
    /// Skip branches whose radius exceeds this. Use to keep leaves off the trunk.
    float maxRadius = 0.05f;
    /// Skip segments whose `depth` is below this.
    int   minDepth = 1;
    /// Only place leaves on segments with no children (chain tips).
    bool  terminalOnly = false;

    /// Average leaves per unit of segment length.
    float perUnitLength  = 20.0f;
    /// 0 = uniform along segment; >0 biases samples toward the tip
    /// (t = 1 - (1-u)^(1+falloff)).
    float densityFalloff = 0.0f;

    /// 0 = leaf forward points purely radially-outward from the branch;
    /// 1 = leaf forward forced toward world up (phototropism).
    float upBias = 0.5f;
    /// ± radians of random pitch around the branch tangent.
    float tiltJitter = 0.3f;
    /// ± radians of random roll around the leaf forward axis.
    float rollJitter = 0.2f;

    /// Base uniform scale applied to every leaf.
    float baseScale = 1.0f;
    /// ± fraction of `baseScale` applied per leaf.
    float scaleJitter = 0.2f;
    /// 0 = ignore branch radius; 1 = scale linearly with (radius / maxRadius).
    /// Lets thinner twigs carry smaller leaves.
    float scaleByRadius = 0.0f;

    /// Minimum distance between accepted leaf origins. 0 disables.
    float dedupRadius = 0.0f;

    /// Deterministic seed.
    uint64_t seed = 0;
};

/// Flat per-leaf instance buffer. `transforms` stride is 16 floats (column-major
/// 4x4 of T * R * uniform-S). `branchRadius` and `branchDepth` are 1 entry per
/// leaf, useful for shading variation downstream.
struct LeafPlacements {
    std::vector<float>   transforms;
    std::vector<float>   branchRadius;
    std::vector<int32_t> branchDepth;
    size_t count() const { return transforms.size() / 16; }
};

/// Compute leaf instance transforms along branch segments.
///
/// Per-leaf local frame matches `leafCard` output: +Z = tip direction,
/// +Y = card normal, +X = side. The returned transform takes a leaf in that
/// local space into world space.
LeafPlacements placeLeavesOnBranches(
    const std::vector<BranchSegment>& segments,
    const LeafPlacementOptions& opts = {});

/// Stamp `leaf` (in its local space) at every placement and return a single
/// merged mesh. Positions are transformed by the full 4x4; normals by the
/// rotation part (uniform scale, so no inverse-transpose needed).
MeshData scatterLeaves(
    const std::vector<BranchSegment>& segments,
    const MeshData& leaf,
    const LeafPlacementOptions& opts = {});

} // namespace bromesh
