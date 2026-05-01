#include "bromesh/procedural/plants/vine.h"

#include "plant_common.h"

#include <algorithm>
#include <cmath>
#include <random>

namespace bromesh {

PlantResult buildVine(const VineParams& params) {
    PlantResult result;
    using namespace plant_internal;

    std::mt19937_64 rng(params.seed);
    std::uniform_real_distribution<float> uni(-1.0f, 1.0f);

    // age01 truncates the path; helix winding is preserved up to the cut.
    const float age01 = std::clamp(params.age01, 0.0f, 1.0f);
    const float effLength = params.length * std::max(0.05f, age01);
    const int samples = std::max(16, static_cast<int>(effLength * 12.0f));
    const float two_pi = 6.28318530717958647692f;
    std::vector<Vec3> path;
    path.reserve(samples + 1);
    for (int i = 0; i <= samples; ++i) {
        float t = static_cast<float>(i) / static_cast<float>(samples);
        float a = t * params.turns * two_pi;
        Vec3 p{
            std::cos(a) * params.helixRadius + uni(rng) * 0.04f,
            t * effLength,
            std::sin(a) * params.helixRadius + uni(rng) * 0.04f
        };
        path.push_back(p);
    }

    auto profile = circleProfile(6, 1.0f);
    SweepOptions opts;
    opts.closeProfile = true;
    opts.capStart = true;
    opts.capEnd = true;
    opts.miterJoints = true;
    opts.profileScale.assign(path.size(), params.radius);
    result.branchMesh = sweep(profile, path, opts);

    // Periodic leaves along the path.
    const int leafCount = std::max(1, static_cast<int>(effLength * params.leafDensity));
    for (int i = 0; i < leafCount; ++i) {
        float t = (static_cast<float>(i) + 0.5f) / static_cast<float>(leafCount);
        size_t idx = static_cast<size_t>(t * (path.size() - 1));
        Vec3 p = path[idx];
        Vec3 fwd = (idx + 1 < path.size()) ? vnormOr(path[idx + 1] - p, {0,1,0})
                                            : Vec3{0, 1, 0};
        LeafInstance L;
        L.position = p;
        L.orientation = quatLookDir(fwd);
        L.scale = 0.12f;
        L.variantIndex = i & 1;
        result.leaves.push_back(L);
    }

    if (!result.branchMesh.empty())
        aabbFromMesh(result.branchMesh, result.aabbMin, result.aabbMax);
    else { result.aabbMin = {0,0,0}; result.aabbMax = {0,effLength,0}; }
    for (const auto& L : result.leaves) {
        Vec3 r{L.scale, L.scale, L.scale};
        updateAabb(result.aabbMin, result.aabbMax, L.position - r);
        updateAabb(result.aabbMin, result.aabbMax, L.position + r);
    }
    return result;
}

} // namespace bromesh
