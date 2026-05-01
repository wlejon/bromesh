#include "bromesh/procedural/plants/succulent.h"

#include "plant_common.h"

#include <algorithm>
#include <cmath>
#include <random>

namespace bromesh {

PlantResult buildSucculent(const SucculentParams& params) {
    PlantResult result;
    using namespace plant_internal;

    const int leaves = std::max(3, params.leafCount);
    std::mt19937_64 rng(params.seed);
    std::uniform_real_distribution<float> uni(0.0f, 1.0f);

    std::vector<MeshData> parts;
    parts.reserve(leaves);

    // Phyllotactic golden angle in radians.
    const float golden = 2.39996323f;

    for (int i = 0; i < leaves; ++i) {
        float a = static_cast<float>(i) * golden;
        // Outer leaves slightly tilted up; inner leaves more vertical.
        float t = static_cast<float>(i) / static_cast<float>(leaves);
        float tilt = params.tilt * (0.4f + 0.7f * t);
        float length = params.leafLength * (0.7f + 0.4f * t + uni(rng) * 0.1f);

        Vec3 outDir{ std::cos(a) * std::cos(tilt), std::sin(tilt), std::sin(a) * std::cos(tilt) };
        outDir = vnorm(outDir);

        std::vector<Vec3> path;
        const int segs = 8;
        for (int s = 0; s <= segs; ++s) {
            float u = static_cast<float>(s) / static_cast<float>(segs);
            // Slight upward curl at the tip.
            Vec3 p = outDir * (u * length);
            p.y += u * u * length * 0.25f;
            path.push_back(p);
        }

        auto profile = diamondProfile(params.leafWidth, params.leafThickness);
        SweepOptions opts;
        opts.closeProfile = true;
        opts.capStart = true;
        opts.capEnd = true;
        opts.miterJoints = true;
        opts.profileScale.resize(path.size());
        for (size_t s = 0; s < path.size(); ++s) {
            float u = static_cast<float>(s) / static_cast<float>(path.size() - 1);
            // Wide near base, taper to a point at the tip.
            float bulge = std::sin(3.14159265f * u);
            opts.profileScale[s] = std::max(0.05f, 0.4f + 0.6f * bulge - 0.5f * u);
        }
        MeshData leaf = sweep(profile, path, opts);
        if (!leaf.empty()) parts.push_back(std::move(leaf));
    }

    if (!parts.empty()) result.branchMesh = mergeMeshes(parts);
    if (!result.branchMesh.empty())
        aabbFromMesh(result.branchMesh, result.aabbMin, result.aabbMax);
    else { result.aabbMin = {0,0,0}; result.aabbMax = {0, params.leafLength, 0}; }
    return result;
}

} // namespace bromesh
