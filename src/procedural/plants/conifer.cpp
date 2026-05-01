#include "bromesh/procedural/plants/conifer.h"

#include "plant_common.h"

#include <algorithm>
#include <cmath>
#include <random>

namespace bromesh {

PlantResult buildConifer(const ConiferParams& params) {
    PlantResult result;
    using namespace plant_internal;

    const float age = std::clamp(params.age01, 0.05f, 1.0f);
    const float H = params.height * age;
    const int whorls = std::max(2, params.whorlCount);
    const int perWhorl = std::max(3, params.branchesPerWhorl);

    std::mt19937_64 rng(params.seed);
    std::uniform_real_distribution<float> uni(0.0f, 1.0f);

    // Trunk path: straight, with a slight noisy lean.
    std::vector<Vec3> trunkPath;
    const int trunkSamples = 16;
    for (int i = 0; i <= trunkSamples; ++i) {
        float t = static_cast<float>(i) / static_cast<float>(trunkSamples);
        Vec3 p{
            (uni(rng) - 0.5f) * 0.04f * H,
            t * H,
            (uni(rng) - 0.5f) * 0.04f * H
        };
        trunkPath.push_back(p);
    }
    auto trunkProfile = circleProfile(10, 1.0f);
    SweepOptions trunkOpts;
    trunkOpts.closeProfile = true;
    trunkOpts.capStart = true;
    trunkOpts.capEnd = true;
    trunkOpts.miterJoints = true;
    trunkOpts.profileScale.resize(trunkPath.size());
    for (size_t i = 0; i < trunkPath.size(); ++i) {
        float t = static_cast<float>(i) / static_cast<float>(trunkPath.size() - 1);
        trunkOpts.profileScale[i] = params.baseRadius * (1.0f - 0.85f * t);
    }
    MeshData trunkMesh = sweep(trunkProfile, trunkPath, trunkOpts);

    // Branches at whorls.
    std::vector<MeshData> parts;
    parts.reserve(whorls * perWhorl + 1);
    parts.push_back(std::move(trunkMesh));

    auto branchProfile = circleProfile(6, 1.0f);
    const float two_pi = 6.28318530717958647692f;
    for (int w = 0; w < whorls; ++w) {
        float tWhorl = (static_cast<float>(w) + 0.5f) / static_cast<float>(whorls);
        float yBase = tWhorl * H * 0.95f + H * 0.05f;
        float trunkR = params.baseRadius * (1.0f - 0.85f * tWhorl);
        // Branches taper from long at the bottom to short at the top.
        float branchLen = (1.0f - tWhorl) * H * 0.55f + 0.2f;
        float angleOffset = uni(rng) * two_pi;
        for (int b = 0; b < perWhorl; ++b) {
            float a = angleOffset + two_pi * static_cast<float>(b) / static_cast<float>(perWhorl);
            Vec3 dir{ std::cos(a), -0.15f, std::sin(a) };
            dir = vnorm(dir);
            Vec3 from{ trunkR * std::cos(a), yBase, trunkR * std::sin(a) };
            // Build a curved branch with a slight droop.
            std::vector<Vec3> path;
            const int segs = 8;
            for (int s = 0; s <= segs; ++s) {
                float t = static_cast<float>(s) / static_cast<float>(segs);
                float droop = -t * t * 0.25f * branchLen;
                Vec3 p = from + dir * (t * branchLen);
                p.y += droop;
                path.push_back(p);
            }
            SweepOptions bopts;
            bopts.closeProfile = true;
            bopts.capStart = false;
            bopts.capEnd = true;
            bopts.miterJoints = true;
            bopts.profileScale.resize(path.size());
            float r0 = std::max(0.012f, trunkR * 0.25f);
            for (size_t s = 0; s < path.size(); ++s) {
                float t = static_cast<float>(s) / static_cast<float>(path.size() - 1);
                bopts.profileScale[s] = r0 * (1.0f - 0.85f * t);
            }
            MeshData bm = sweep(branchProfile, path, bopts);
            if (!bm.empty()) parts.push_back(std::move(bm));

            // Foliage: leaves along the back half of the branch.
            int leafN = std::max(4, static_cast<int>(branchLen * 14.0f));
            for (int i = 0; i < leafN; ++i) {
                float t = 0.4f + 0.6f * static_cast<float>(i) / static_cast<float>(leafN - 1);
                LeafInstance L;
                Vec3 p = from + dir * (t * branchLen);
                p.y += -t * t * 0.25f * branchLen;
                p.x += (uni(rng) - 0.5f) * 0.05f;
                p.z += (uni(rng) - 0.5f) * 0.05f;
                L.position = p;
                L.orientation = quatLookDir(dir);
                L.scale = 0.05f * (1.0f - t * 0.5f);
                L.variantIndex = 0;
                result.leaves.push_back(L);
            }
        }
    }

    result.branchMesh = mergeMeshes(parts);
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
