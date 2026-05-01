#include "bromesh/procedural/plants/fern.h"

#include "bromesh/procedural/lsystem.h"
#include "plant_common.h"

#include <algorithm>
#include <cmath>

namespace bromesh {

PlantResult buildFern(const FernParams& params) {
    PlantResult result;
    using namespace plant_internal;

    const int pairs = std::max(2, params.leafletPairs);
    // Same seed -> same mature frond. age01 truncates the rachis from the
    // tip by scaling total length; the L-system structure is fixed.
    const float age01 = std::clamp(params.age01, 0.05f, 1.0f);

    // L-system: A is a frond unit. After N iterations the rachis is a
    // chain of F segments interspersed with paired leaflet markers L+ / L-.
    LSystem ls;
    ls.setAxiom({ {'A', {1.0f}} });

    ProductionRule grow;
    grow.predecessor = 'A';
    grow.weight = 1.0f;
    grow.successor = [](const std::vector<float>& p) {
        float t = p.empty() ? 1.0f : p[0];
        return std::vector<Module>{
            {'F', {t}},
            {'[', {}},
            {'+', {1.0f}}, {'L', {t}},
            {']', {}},
            {'[', {}},
            {'-', {1.0f}}, {'L', {t}},
            {']', {}},
            {'A', { t * 0.95f }}
        };
    };
    ls.addRule(grow);

    auto modules = ls.derive(pairs, params.seed);

    // Walk modules to lay down a curved rachis path and place leaflets.
    Vec3 cur{0, 0, 0};
    Vec3 dir{0, 1, 0};
    const float total = params.length * age01;
    // Step length so the F count along the rachis spans `total`.
    int fCount = 0;
    for (const auto& m : modules) if (m.symbol == 'F') ++fCount;
    if (fCount == 0) fCount = 1;
    const float step = total / static_cast<float>(fCount);
    const float bendPerStep = params.curvature / static_cast<float>(fCount);

    std::vector<Vec3> rachis;
    rachis.push_back(cur);

    std::vector<MeshData> parts;
    parts.reserve(modules.size());

    int leafletIndex = 0;
    int side = 1; // +1 / -1 alternated by '+' / '-' markers.
    Vec3 sideDir{1, 0, 0};

    for (const auto& m : modules) {
        if (m.symbol == 'F') {
            // Rotate dir slightly toward forward-down (curl).
            float c = std::cos(bendPerStep);
            float s = std::sin(bendPerStep);
            // Rotate around X axis: y' = c*y - s*z; z' = s*y + c*z. Path
            // bends along world-Z so the frond curls forward.
            Vec3 nd{ dir.x, c * dir.y - s * dir.z, s * dir.y + c * dir.z };
            dir = vnorm(nd);
            cur = cur + dir * step;
            rachis.push_back(cur);
        } else if (m.symbol == '+') {
            side = +1;
        } else if (m.symbol == '-') {
            side = -1;
        } else if (m.symbol == 'L') {
            // Leaflet at current position, perpendicular to dir.
            float t = static_cast<float>(leafletIndex) / static_cast<float>(std::max(1, fCount));
            float taper = std::sin(3.14159265f * t);
            float length = params.leafletLength * std::max(0.1f, taper);
            Vec3 ortho = vnorm(vcross(dir, {0, 0, 1}));
            if (vdot(ortho, ortho) < 0.5f) ortho = vnorm(vcross(dir, {1, 0, 0}));
            sideDir = ortho * static_cast<float>(side);

            std::vector<Vec3> path;
            const int ls_segs = 6;
            for (int i = 0; i <= ls_segs; ++i) {
                float u = static_cast<float>(i) / static_cast<float>(ls_segs);
                Vec3 p = cur + sideDir * (u * length) + dir * (u * length * 0.1f);
                path.push_back(p);
            }
            auto profile = diamondProfile(0.025f * std::max(0.2f, taper),
                                          0.005f);
            SweepOptions opts;
            opts.closeProfile = true;
            opts.capStart = false;
            opts.capEnd = true;
            opts.miterJoints = false;
            opts.profileScale.resize(path.size());
            for (size_t i = 0; i < path.size(); ++i) {
                float u = static_cast<float>(i) / static_cast<float>(path.size() - 1);
                opts.profileScale[i] = std::max(0.1f, 1.0f - u);
            }
            MeshData lm = sweep(profile, path, opts);
            if (!lm.empty()) parts.push_back(std::move(lm));
            ++leafletIndex;
        }
        // '[' and ']' ignored — leaflets are absolute placements.
    }

    // Mesh the rachis as a tapered tube.
    if (rachis.size() >= 2) {
        auto profile = circleProfile(6, 1.0f);
        SweepOptions opts;
        opts.closeProfile = true;
        opts.capStart = true;
        opts.capEnd = true;
        opts.miterJoints = true;
        opts.profileScale.resize(rachis.size());
        for (size_t i = 0; i < rachis.size(); ++i) {
            float t = static_cast<float>(i) / static_cast<float>(rachis.size() - 1);
            opts.profileScale[i] = params.stemRadius * std::max(0.15f, 1.0f - t);
        }
        MeshData stem = sweep(profile, rachis, opts);
        if (!stem.empty()) parts.push_back(std::move(stem));
    }

    if (!parts.empty()) result.branchMesh = mergeMeshes(parts);
    if (!result.branchMesh.empty())
        aabbFromMesh(result.branchMesh, result.aabbMin, result.aabbMax);
    else { result.aabbMin = {0,0,0}; result.aabbMax = {0, params.length, 0}; }
    return result;
}

} // namespace bromesh
