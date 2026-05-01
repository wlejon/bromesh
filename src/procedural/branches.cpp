#include "bromesh/procedural/branches.h"

#include "bromesh/manipulation/merge.h"
#include "bromesh/manipulation/normals.h"
#include "bromesh/manipulation/sweep.h"
#include "bromesh/procedural/vec_math.h"

#include <cmath>

namespace bromesh {

static std::vector<Vec2> circleProfile(int segments, float radius) {
    std::vector<Vec2> out;
    out.reserve(segments);
    const float two_pi = 6.28318530717958647692f;
    for (int i = 0; i < segments; ++i) {
        float a = two_pi * static_cast<float>(i) / static_cast<float>(segments);
        out.push_back({std::cos(a) * radius, std::sin(a) * radius});
    }
    return out;
}

MeshData meshBranches(const std::vector<BranchSegment>& segs, int sides) {
    if (segs.empty()) return {};

    std::vector<std::vector<int>> children(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) {
        int p = segs[i].parent;
        if (p >= 0 && static_cast<size_t>(p) < segs.size()) {
            children[p].push_back(static_cast<int>(i));
        }
    }

    auto isChainStart = [&](int i) -> bool {
        int p = segs[i].parent;
        if (p < 0) return true;
        return children[p].size() != 1;
    };

    auto profile = circleProfile(sides, 1.0f);
    std::vector<MeshData> parts;
    parts.reserve(segs.size());

    for (size_t i = 0; i < segs.size(); ++i) {
        if (!isChainStart(static_cast<int>(i))) continue;

        std::vector<Vec3> path;
        std::vector<float> radii;

        int p0 = segs[i].parent;
        float r0 = (p0 >= 0 && segs[p0].radius > 0.0f) ? segs[p0].radius
                  : (segs[i].radius > 0.0f ? segs[i].radius : 0.02f);
        path.push_back(segs[i].from);
        radii.push_back(r0);

        int cur = static_cast<int>(i);
        while (true) {
            const BranchSegment& s = segs[cur];
            if (vdist2(s.from, s.to) > 1e-12f) {
                path.push_back(s.to);
                radii.push_back(s.radius > 0.0f ? s.radius : radii.back());
            }
            if (children[cur].size() != 1) break;
            cur = children[cur][0];
        }

        if (path.size() < 2) continue;

        SweepOptions opts;
        opts.closeProfile = true;
        opts.capStart = (segs[i].parent < 0);
        opts.capEnd   = children[cur].empty();
        opts.miterJoints = false;
        opts.profileScale = radii;

        MeshData part = sweep(profile, path, opts);
        if (!part.empty()) parts.push_back(std::move(part));
    }

    if (parts.empty()) return {};
    MeshData merged = mergeMeshes(parts);
    computeNormals(merged);
    return merged;
}

} // namespace bromesh
