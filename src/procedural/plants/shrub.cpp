#include "bromesh/procedural/plants/shrub.h"

#include "plant_common.h"

#include <algorithm>
#include <cmath>
#include <random>

namespace bromesh {

PlantResult buildShrub(const ShrubParams& params) {
    PlantResult result;
    using namespace plant_internal;

    // Mature structure is fixed; age01 trims by depth so the same seed
    // produces a consistent skeleton growing outward.
    const float age01 = std::clamp(params.age01, 0.0f, 1.0f);
    const float H = params.height;
    const float R = params.radius;
    const int stems = std::max(1, params.stemCount);
    const int attractorCount = std::max(8, params.attractorCount);

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

    // Filter by depth-based birth time. Same seed -> same mature skeleton
    // for all ages; age01 just trims descendants beyond the current front.
    int maxDepth = 1;
    for (const auto& s : segs) if (s.depth > maxDepth) maxDepth = s.depth;
    auto h01 = [&](int idx)->float {
        uint64_t x = params.seed ^ (uint64_t)(uint32_t)idx * 0x9E3779B97F4A7C15ull;
        x ^= x >> 33; x *= 0xff51afd7ed558ccdull;
        x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ull;
        x ^= x >> 33;
        return static_cast<float>((x >> 40) & 0xFFFFFFu) / static_cast<float>(0x1000000u);
    };
    std::vector<bool> alive(segs.size(), true);
    if (age01 < 1.0f) {
        for (size_t i = 0; i < segs.size(); ++i) {
            float dt = static_cast<float>(segs[i].depth) / static_cast<float>(maxDepth);
            float jitter = (h01(static_cast<int>(i)) - 0.5f) * 0.05f;
            alive[i] = std::clamp(dt + jitter, 0.0f, 1.0f) <= age01;
        }
        for (size_t i = 0; i < segs.size(); ++i) {
            int p = segs[i].parent;
            if (p >= 0 && !alive[p]) alive[i] = false;
        }
    }
    std::vector<int> remap(segs.size(), -1);
    std::vector<BranchSegment> kept;
    kept.reserve(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) {
        if (!alive[i]) continue;
        BranchSegment s = segs[i];
        s.parent = (s.parent >= 0) ? remap[s.parent] : -1;
        remap[i] = static_cast<int>(kept.size());
        kept.push_back(s);
    }

    result.branchMesh = meshBranches(kept, 6);

    // Leaves on terminals of the surviving structure.
    std::vector<int> childCount(kept.size(), 0);
    for (const auto& s : kept) if (s.parent >= 0) ++childCount[s.parent];
    for (size_t i = 0; i < kept.size(); ++i) {
        if (childCount[i] != 0) continue;
        if (kept[i].parent == -1) continue;
        LeafInstance L;
        L.position = kept[i].to;
        Vec3 fwd = vnormOr(kept[i].to - kept[i].from, {0, 1, 0});
        L.orientation = quatLookDir(fwd);
        L.scale = 0.07f * (0.7f + h01(static_cast<int>(i) ^ 0x55) * 0.6f);
        L.variantIndex = static_cast<int>(
            (uint32_t)(h01(static_cast<int>(i) ^ 0xAA) * 4.0f)) & 3;
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
