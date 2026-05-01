#include "bromesh/procedural/plants/tree.h"

#include "plant_common.h"

#include <cmath>
#include <random>

namespace bromesh {

PlantResult buildTree(const TreeParams& params) {
    PlantResult result;
    using namespace plant_internal;

    const float age = std::clamp(params.age01, 0.05f, 1.0f);
    const float H = params.height * age;
    const float CR = params.canopyRadius * age;
    const int attractorCount = std::max(8, static_cast<int>(params.attractorCount * age));

    // Generate hemispheric attractor cloud above the trunk.
    std::mt19937_64 rng(params.seed);
    std::uniform_real_distribution<float> uni(0.0f, 1.0f);
    const float canopyBase = H * 0.45f;
    const float canopyTop = H * 1.05f;
    std::vector<Vec3> attractors;
    attractors.reserve(attractorCount);
    for (int i = 0; i < attractorCount; ++i) {
        // Rejection-sample a point in a sphere then squish vertically.
        Vec3 p;
        do {
            p = {
                uni(rng) * 2.0f - 1.0f,
                uni(rng) * 2.0f - 1.0f,
                uni(rng) * 2.0f - 1.0f
            };
        } while (vdot(p, p) > 1.0f);
        Vec3 a = {
            p.x * CR,
            canopyBase + (uni(rng) * 0.5f + 0.5f) * (canopyTop - canopyBase),
            p.z * CR
        };
        attractors.push_back(a);
    }

    SpaceColonizationOptions opts;
    opts.attractionRadius = CR * 0.7f;
    opts.killRadius = std::max(params.trunkRadius * 1.5f, CR * 0.08f);
    opts.segmentLength = std::max(0.1f, H * 0.05f);
    opts.maxIterations = 220;
    opts.tropism = {0, 1, 0};
    opts.tropismWeight = 0.35f;

    // Pre-grow a straight trunk up to the canopy base so colonization has
    // a starting node within reach of the attractor cloud. The seed point
    // for colonization is the top of this trunk; the trunk segments are
    // inserted up front so the segment list still forms a tree.
    std::vector<BranchSegment> segs;
    {
        const float trunkTop = canopyBase * 0.9f;
        const float step = opts.segmentLength;
        int trunkSegs = std::max(1, static_cast<int>(trunkTop / step));
        Vec3 prev{0, 0, 0};
        // Root marker.
        BranchSegment root;
        root.parent = -1;
        root.from = prev;
        root.to = prev;
        root.depth = 0;
        segs.push_back(root);
        for (int i = 1; i <= trunkSegs; ++i) {
            BranchSegment s;
            s.parent = static_cast<int>(segs.size()) - 1;
            s.from = prev;
            s.to = { 0.0f, static_cast<float>(i) * step, 0.0f };
            s.depth = i;
            segs.push_back(s);
            prev = s.to;
        }
    }

    // Continue growth from the top of the trunk using space colonization.
    // We re-run the algorithm by feeding it the existing tree as seeds and
    // splicing the result in. The simplest path: just call spaceColonize
    // with the trunk top as the seed and concatenate, fixing parent indices.
    {
        Vec3 trunkTopPos = segs.back().to;
        int trunkTopSeg = static_cast<int>(segs.size()) - 1;
        std::vector<BranchSegment> grown = spaceColonize(
            attractors, { trunkTopPos }, {0, 1, 0}, opts);
        // grown[0] is the synthetic root segment for trunkTopPos with parent
        // == -1 and from == to. Skip it; reparent its direct children to
        // trunkTopSeg.
        const int offset = static_cast<int>(segs.size()) - 1; // we drop grown[0]
        for (size_t i = 1; i < grown.size(); ++i) {
            BranchSegment s = grown[i];
            if (s.parent == 0) s.parent = trunkTopSeg;
            else if (s.parent > 0) s.parent = s.parent + offset;
            s.depth += static_cast<int>(segs.size()) - 1;
            segs.push_back(s);
        }
    }

    thickenBranches(segs, std::max(0.01f, params.trunkRadius * 0.12f), 2.5f);

    // Boost root segment radius to user-requested trunk radius.
    for (auto& s : segs) {
        if (s.parent == -1) s.radius = std::max(s.radius, params.trunkRadius);
    }

    result.branchMesh = meshBranches(segs, 8);

    // Place leaves on terminal segments (leaves: segments with no children).
    std::vector<int> childCount(segs.size(), 0);
    for (const auto& s : segs) if (s.parent >= 0) ++childCount[s.parent];
    const float leafScaleBase = std::max(0.05f, H * 0.04f);
    for (size_t i = 0; i < segs.size(); ++i) {
        if (childCount[i] != 0) continue;
        if (segs[i].parent == -1) continue;
        LeafInstance L;
        L.position = segs[i].to;
        Vec3 fwd = vnormOr(segs[i].to - segs[i].from, {0, 1, 0});
        L.orientation = quatLookDir(fwd);
        L.scale = leafScaleBase * (0.7f + uni(rng) * 0.6f);
        L.variantIndex = static_cast<int>(rng() & 3);
        result.leaves.push_back(L);
    }

    if (!result.branchMesh.empty()) {
        plant_internal::aabbFromMesh(result.branchMesh, result.aabbMin, result.aabbMax);
    } else {
        result.aabbMin = {0, 0, 0};
        result.aabbMax = {0, H, 0};
    }
    for (const auto& L : result.leaves) {
        Vec3 r{L.scale, L.scale, L.scale};
        plant_internal::updateAabb(result.aabbMin, result.aabbMax, L.position - r);
        plant_internal::updateAabb(result.aabbMin, result.aabbMax, L.position + r);
    }
    return result;
}

} // namespace bromesh
