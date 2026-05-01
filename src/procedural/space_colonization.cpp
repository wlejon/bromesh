#include "bromesh/procedural/space_colonization.h"

#include "bromesh/optimization/spatial_hash.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <vector>

namespace bromesh {

std::vector<BranchSegment> spaceColonize(
    const std::vector<Vec3>& attractors,
    const std::vector<Vec3>& seedPoints,
    const Vec3& initialDirection,
    const SpaceColonizationOptions& opts) {
    std::vector<BranchSegment> segments;
    if (seedPoints.empty()) return segments;

    // Tree nodes: position + parent segment index (-1 for seed roots).
    struct Node {
        Vec3 pos;
        int seg;   // index into `segments` of the segment ending at this node
        int depth;
    };
    std::vector<Node> nodes;
    nodes.reserve(seedPoints.size() + 64);

    // Spatial hash of node positions.
    float cell = std::max(opts.attractionRadius, opts.segmentLength) * 1.0f;
    if (cell <= 0.0f) cell = 1.0f;
    SpatialHash3D nodeHash(cell);

    Vec3 initDir = vnormOr(initialDirection, {0, 1, 0});
    for (const Vec3& sp : seedPoints) {
        BranchSegment root;
        root.parent = -1;
        root.from = sp;
        root.to = sp;
        root.depth = 0;
        segments.push_back(root);
        Node n{sp, static_cast<int>(segments.size()) - 1, 0};
        nodeHash.insert(sp, static_cast<int32_t>(nodes.size()));
        nodes.push_back(n);
    }

    // Live attractors: index list into `attractors`.
    std::vector<bool> alive(attractors.size(), true);

    const float r2Attr = opts.attractionRadius * opts.attractionRadius;
    const float kill2 = opts.killRadius * opts.killRadius;

    std::vector<Vec3> growthAccum(0);
    std::vector<int>  growthCount(0);

    for (int iter = 0; iter < opts.maxIterations; ++iter) {
        growthAccum.assign(nodes.size(), Vec3{0, 0, 0});
        growthCount.assign(nodes.size(), 0);

        // For each live attractor, find nearest node within attractionRadius.
        bool anyInfluence = false;
        for (size_t a = 0; a < attractors.size(); ++a) {
            if (!alive[a]) continue;
            int32_t nearest = nodeHash.nearest(attractors[a], opts.attractionRadius);
            if (nearest < 0) continue;
            (void)r2Attr;
            Vec3 dir = vnorm(attractors[a] - nodes[nearest].pos);
            growthAccum[nearest] += dir;
            growthCount[nearest] += 1;
            anyInfluence = true;
        }
        if (!anyInfluence) break;

        // Grow new nodes from any node that received attractor pull.
        size_t newStart = nodes.size();
        for (size_t i = 0; i < newStart; ++i) {
            if (growthCount[i] == 0) continue;
            Vec3 dir = vnorm(growthAccum[i]);
            if (opts.tropismWeight > 0.0f) {
                dir = vnormOr(dir + opts.tropism * opts.tropismWeight, dir);
            }
            (void)initDir;
            Vec3 newPos = nodes[i].pos + dir * opts.segmentLength;

            BranchSegment seg;
            // The segment ending at the parent node also identifies the
            // parent in the segment list — except for the synthetic root
            // segments where from == to. Use that segment as the parent so
            // the segment list forms a tree.
            seg.parent = nodes[i].seg;
            seg.from = nodes[i].pos;
            seg.to = newPos;
            seg.depth = nodes[i].depth + 1;
            segments.push_back(seg);

            Node nn{newPos, static_cast<int>(segments.size()) - 1, nodes[i].depth + 1};
            nodeHash.insert(newPos, static_cast<int32_t>(nodes.size()));
            nodes.push_back(nn);
        }

        // Kill attractors near any node (including the brand-new ones).
        for (size_t a = 0; a < attractors.size(); ++a) {
            if (!alive[a]) continue;
            int32_t hit = nodeHash.nearest(attractors[a], opts.killRadius);
            if (hit >= 0) {
                Vec3 d = attractors[a] - nodes[hit].pos;
                if (vdot(d, d) <= kill2) alive[a] = false;
            }
        }

        if (newStart == nodes.size()) break;
    }

    return segments;
}

void thickenBranches(std::vector<BranchSegment>& segments,
                     float leafRadius,
                     float pipeExp) {
    if (segments.empty()) return;
    if (pipeExp <= 0.0f) pipeExp = 2.5f;

    // Build child lists.
    std::vector<std::vector<int>> children(segments.size());
    for (size_t i = 0; i < segments.size(); ++i) {
        int p = segments[i].parent;
        if (p >= 0 && static_cast<size_t>(p) < segments.size()) {
            children[p].push_back(static_cast<int>(i));
        }
    }
    // Post-order traversal: process children before parents. Iterative.
    std::vector<int> order;
    order.reserve(segments.size());
    std::vector<int> stack;
    std::vector<bool> visited(segments.size(), false);
    for (size_t i = 0; i < segments.size(); ++i) {
        if (segments[i].parent == -1) stack.push_back(static_cast<int>(i));
    }
    while (!stack.empty()) {
        int idx = stack.back();
        if (!visited[idx]) {
            visited[idx] = true;
            for (int c : children[idx]) stack.push_back(c);
        } else {
            stack.pop_back();
            // Already on order? No — we use a different scheme.
        }
    }
    // Simpler: build pre-order then reverse for post-order processing.
    order.clear();
    stack.clear();
    std::fill(visited.begin(), visited.end(), false);
    for (size_t i = 0; i < segments.size(); ++i) {
        if (segments[i].parent == -1) stack.push_back(static_cast<int>(i));
    }
    while (!stack.empty()) {
        int idx = stack.back();
        stack.pop_back();
        if (visited[idx]) continue;
        visited[idx] = true;
        order.push_back(idx);
        for (int c : children[idx]) stack.push_back(c);
    }
    std::reverse(order.begin(), order.end());

    for (int idx : order) {
        if (children[idx].empty()) {
            segments[idx].radius = leafRadius;
        } else {
            float sum = 0.0f;
            for (int c : children[idx]) {
                sum += std::pow(segments[c].radius, pipeExp);
            }
            segments[idx].radius = std::pow(sum, 1.0f / pipeExp);
        }
    }
}

} // namespace bromesh
