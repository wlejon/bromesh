#include "bromesh/procedural/plants/shrub.h"

#include "plant_common.h"

#include <algorithm>
#include <cmath>
#include <random>

namespace bromesh {

PlantResult buildShrub(const ShrubParams& params) {
    PlantResult result;
    using namespace plant_internal;

    const float age = std::clamp(params.age01, 0.05f, 1.0f);
    const float H = params.height * age;
    const float R = params.radius * age;
    const int stems = std::max(1, params.stemCount);
    const int attractorCount = std::max(8, static_cast<int>(params.attractorCount * age));

    std::mt19937_64 rng(params.seed);
    std::uniform_real_distribution<float> uni(0.0f, 1.0f);

    // Cluster of seed points at the base.
    std::vector<Vec3> seeds;
    seeds.reserve(stems);
    const float two_pi = 6.28318530717958647692f;
    for (int i = 0; i < stems; ++i) {
        float a = two_pi * static_cast<float>(i) / static_cast<float>(stems);
        float r = R * 0.1f * (uni(rng) + 0.2f);
        seeds.push_back({ std::cos(a) * r, 0.0f, std::sin(a) * r });
    }

    // Attractor cloud: dense and low.
    std::vector<Vec3> attractors;
    attractors.reserve(attractorCount);
    for (int i = 0; i < attractorCount; ++i) {
        Vec3 p;
        do {
            p = { uni(rng)*2-1, uni(rng)*2-1, uni(rng)*2-1 };
        } while (vdot(p, p) > 1.0f);
        attractors.push_back({ p.x * R, H * (0.2f + 0.8f * uni(rng)), p.z * R });
    }

    SpaceColonizationOptions opts;
    opts.attractionRadius = R * 0.6f;
    opts.killRadius = std::max(0.05f, R * 0.08f);
    opts.segmentLength = std::max(0.06f, H * 0.07f);
    opts.maxIterations = 160;
    opts.tropism = {0, 1, 0};
    opts.tropismWeight = 0.25f;

    auto segs = spaceColonize(attractors, seeds, {0, 1, 0}, opts);
    thickenBranches(segs, 0.012f, 2.5f);

    result.branchMesh = meshBranches(segs, 6);

    // Leaves: terminals + a fraction of interior tips.
    std::vector<int> childCount(segs.size(), 0);
    for (const auto& s : segs) if (s.parent >= 0) ++childCount[s.parent];
    for (size_t i = 0; i < segs.size(); ++i) {
        if (childCount[i] != 0) continue;
        if (segs[i].parent == -1) continue;
        LeafInstance L;
        L.position = segs[i].to;
        Vec3 fwd = vnormOr(segs[i].to - segs[i].from, {0, 1, 0});
        L.orientation = quatLookDir(fwd);
        L.scale = 0.07f * (0.7f + uni(rng) * 0.6f);
        L.variantIndex = static_cast<int>(rng() & 3);
        result.leaves.push_back(L);
    }

    if (!result.branchMesh.empty())
        aabbFromMesh(result.branchMesh, result.aabbMin, result.aabbMax);
    else { result.aabbMin = {0,0,0}; result.aabbMax = {0,H,0}; }
    for (const auto& L : result.leaves) {
        Vec3 r{L.scale, L.scale, L.scale};
        updateAabb(result.aabbMin, result.aabbMax, L.position - r);
        updateAabb(result.aabbMin, result.aabbMax, L.position + r);
    }
    return result;
}

} // namespace bromesh
