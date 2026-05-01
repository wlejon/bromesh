#pragma once

#include "bromesh/procedural/vec_math.h"

#include <vector>

namespace bromesh {

struct SpaceColonizationOptions {
    /// Radius within which an attractor influences the nearest tree node.
    float attractionRadius = 5.0f;
    /// Attractors closer than this to any node are consumed.
    float killRadius = 0.5f;
    /// Step size per growth iteration.
    float segmentLength = 0.3f;
    /// Maximum growth iterations.
    int maxIterations = 200;
    /// Direction biased into every growth step (e.g. up-vector). Combined
    /// with the attractor-derived direction, weighted by `tropismWeight`.
    Vec3 tropism = {0, 0, 0};
    float tropismWeight = 0.0f;
};

/// One segment of a generated branch tree. Segments index their parent;
/// `parent == -1` marks a root (one of the seed points).
struct BranchSegment {
    int parent = -1;
    Vec3 from{};
    Vec3 to{};
    float radius = 0.0f;
    int depth = 0;
};

/// Run Runions-style space colonization. Returns segments in creation order;
/// each segment's parent index refers to an earlier entry. Roots
/// (one per seed point) have parent = -1 and `from == to == seedPoint`.
std::vector<BranchSegment> spaceColonize(
    const std::vector<Vec3>& attractors,
    const std::vector<Vec3>& seedPoints,
    const Vec3& initialDirection,
    const SpaceColonizationOptions& opts);

/// Compute branch radii using the pipe model. Leaves get `leafRadius`;
/// each parent's radius is the pipe-model sum of its children:
///   r_parent = (sum of r_child^pipeExp)^(1/pipeExp).
/// Modifies segments in place.
void thickenBranches(std::vector<BranchSegment>& segments,
                     float leafRadius = 0.02f,
                     float pipeExp = 2.5f);

} // namespace bromesh
