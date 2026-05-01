#include "bromesh/procedural/plants/grass_tuft.h"

#include "plant_common.h"

#include <algorithm>
#include <cmath>
#include <random>

namespace bromesh {

PlantResult buildGrassTuft(const GrassTuftParams& params) {
    PlantResult result;
    using namespace plant_internal;

    const int blades = std::max(1, params.bladeCount);
    const float age01 = std::clamp(params.age01, 0.05f, 1.0f);
    std::mt19937_64 rng(params.seed);
    std::uniform_real_distribution<float> uni(0.0f, 1.0f);

    std::vector<MeshData> parts;
    parts.reserve(blades);
    const float two_pi = 6.28318530717958647692f;

    for (int i = 0; i < blades; ++i) {
        float ang = two_pi * static_cast<float>(i) / static_cast<float>(blades) + uni(rng) * 0.4f;
        float br = params.baseRadius * (0.4f + 0.6f * uni(rng));
        Vec3 base{ std::cos(ang) * br, 0.0f, std::sin(ang) * br };

        float bladeH = params.height * age01 * (0.7f + 0.6f * uni(rng));
        float tipBend = params.bend * (0.6f + 0.8f * uni(rng));
        Vec3 outDir{ std::cos(ang), 0.0f, std::sin(ang) };

        // Curved arc path.
        std::vector<Vec3> path;
        const int segs = 8;
        for (int s = 0; s <= segs; ++s) {
            float t = static_cast<float>(s) / static_cast<float>(segs);
            float lateral = std::sin(t * tipBend) * bladeH * 0.35f;
            float vertical = std::cos(t * tipBend) * bladeH * t;
            path.push_back({
                base.x + outDir.x * lateral,
                base.y + vertical,
                base.z + outDir.z * lateral
            });
        }

        auto profile = diamondProfile(params.bladeWidth, params.bladeWidth * 0.15f);
        SweepOptions opts;
        opts.closeProfile = true;
        opts.capStart = false;
        opts.capEnd = true;
        opts.miterJoints = true;
        opts.profileScale.resize(path.size());
        for (size_t s = 0; s < path.size(); ++s) {
            float t = static_cast<float>(s) / static_cast<float>(path.size() - 1);
            opts.profileScale[s] = std::max(0.05f, 1.0f - t);
        }
        // Twist the blade so the broad side faces outward initially.
        opts.twist.assign(path.size(), -ang);
        MeshData blade = sweep(profile, path, opts);
        if (!blade.empty()) parts.push_back(std::move(blade));
    }

    if (!parts.empty()) result.branchMesh = mergeMeshes(parts);
    if (!result.branchMesh.empty())
        aabbFromMesh(result.branchMesh, result.aabbMin, result.aabbMax);
    else { result.aabbMin = {0,0,0}; result.aabbMax = {0,params.height,0}; }
    return result;
}

} // namespace bromesh
