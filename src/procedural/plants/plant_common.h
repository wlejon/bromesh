#pragma once

// Internal helpers shared by plant builders. Not part of the public API.

#include "bromesh/manipulation/merge.h"
#include "bromesh/manipulation/normals.h"
#include "bromesh/manipulation/sweep.h"
#include "bromesh/mesh_data.h"
#include "bromesh/procedural/plants/plant_result.h"
#include "bromesh/procedural/space_colonization.h"
#include "bromesh/procedural/vec_math.h"

#include <cmath>
#include <random>
#include <vector>

namespace bromesh {
namespace plant_internal {

inline std::vector<Vec2> circleProfile(int segments, float radius) {
    std::vector<Vec2> out;
    out.reserve(segments);
    const float two_pi = 6.28318530717958647692f;
    for (int i = 0; i < segments; ++i) {
        float a = two_pi * static_cast<float>(i) / static_cast<float>(segments);
        out.push_back({std::cos(a) * radius, std::sin(a) * radius});
    }
    return out;
}

inline std::vector<Vec2> diamondProfile(float halfWidth, float halfThickness) {
    return {
        { halfWidth, 0.0f },
        { 0.0f,  halfThickness },
        { -halfWidth, 0.0f },
        { 0.0f, -halfThickness }
    };
}

/// Mesh a single tapered branch segment as a 2-ring sweep.
inline MeshData meshBranchSweep(Vec3 a, Vec3 b, float ra, float rb, int sides,
                                bool capStart, bool capEnd) {
    MeshData m;
    if (vdist2(a, b) < 1e-12f) return m;
    auto profile = circleProfile(sides, 1.0f);
    SweepOptions opts;
    opts.closeProfile = true;
    opts.capStart = capStart;
    opts.capEnd = capEnd;
    opts.miterJoints = false;
    opts.profileScale = { ra, rb };
    return sweep(profile, { a, b }, opts);
}

/// Mesh an entire branch tree as one merged mesh. Adjacent segments share
/// joint radii (child's start radius == parent's end radius == parent's
/// stored radius), so siblings meet exactly at their common parent. Terminal
/// segments are capped; interior joints are left open since their successors
/// fill the gap.
inline MeshData meshBranches(const std::vector<BranchSegment>& segs, int sides) {
    std::vector<int> childCount(segs.size(), 0);
    for (const auto& s : segs) {
        if (s.parent >= 0 && static_cast<size_t>(s.parent) < childCount.size())
            ++childCount[s.parent];
    }
    std::vector<MeshData> parts;
    parts.reserve(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) {
        const BranchSegment& s = segs[i];
        if (vdist2(s.from, s.to) < 1e-12f) continue;
        float pr = (s.parent >= 0 && static_cast<size_t>(s.parent) < segs.size())
                 ? segs[s.parent].radius : 0.0f;
        float ra = pr > 0.0f ? pr : (s.radius > 0.0f ? s.radius : 0.02f);
        float rb = s.radius > 0.0f ? s.radius : ra;
        bool terminal = (childCount[i] == 0);
        MeshData part = meshBranchSweep(s.from, s.to, ra, rb, sides,
                                        /*capStart*/ false,
                                        /*capEnd*/ terminal);
        if (!part.empty()) parts.push_back(std::move(part));
    }
    if (parts.empty()) return {};
    MeshData merged = mergeMeshes(parts);
    bromesh::computeNormals(merged);
    return merged;
}

inline void updateAabb(Vec3& mn, Vec3& mx, Vec3 p) {
    if (p.x < mn.x) mn.x = p.x; if (p.y < mn.y) mn.y = p.y; if (p.z < mn.z) mn.z = p.z;
    if (p.x > mx.x) mx.x = p.x; if (p.y > mx.y) mx.y = p.y; if (p.z > mx.z) mx.z = p.z;
}

inline void aabbFromMesh(const MeshData& m, Vec3& mn, Vec3& mx) {
    if (m.empty()) {
        mn = {0,0,0}; mx = {0,0,0};
        return;
    }
    mn = { m.positions[0], m.positions[1], m.positions[2] };
    mx = mn;
    const size_t vc = m.vertexCount();
    for (size_t v = 1; v < vc; ++v) {
        Vec3 p{ m.positions[v*3], m.positions[v*3+1], m.positions[v*3+2] };
        updateAabb(mn, mx, p);
    }
}

inline Quat quatLookDir(Vec3 forward) {
    Vec3 f = vnormOr(forward, {0, 1, 0});
    return quatFromTo({0, 1, 0}, f);
}

} // namespace plant_internal
} // namespace bromesh
