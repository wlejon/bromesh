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

/// Mesh an entire branch tree as one merged mesh. Each maximal chain of
/// single-child segments is meshed as one multi-ring sweep so consecutive
/// rings share parallel-transport orientation and vertex positions — no
/// per-segment ring discontinuities. New chains begin at roots and forks;
/// chain endpoints are capped only when terminal (leaf), so forks blend via
/// the parent's continuing radius rather than a flat disc.
inline MeshData meshBranches(const std::vector<BranchSegment>& segs, int sides) {
    if (segs.empty()) return {};

    std::vector<std::vector<int>> children(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) {
        int p = segs[i].parent;
        if (p >= 0 && static_cast<size_t>(p) < segs.size()) {
            children[p].push_back(static_cast<int>(i));
        }
    }

    // A segment starts a new chain iff its parent is a root or a fork (i.e.
    // does not have exactly one child). Otherwise it continues the parent's
    // chain.
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

        // Initial ring sits at the chain's start point with the parent's
        // stored end radius (so a child cylinder begins exactly where the
        // parent's continuing radius leaves off). Falls back to the chain's
        // own first radius for true roots.
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

/// Build a leaf orientation from the branch direction at the leaf's
/// attachment. The leaf mesh is assumed to lie in the local XZ plane with
/// its surface normal as local +Y (the plane primitive convention). We want:
///   - leaf surface facing up-and-outward (catching sun), not edge-on
///   - leaf "stem" axis (local +Z) aligned with the branch tangent
/// so the resulting normal is the bisector of branch-direction and world-up.
/// `quatFromTo({0,1,0}, fwd)` (the prior behaviour) instead pointed the
/// surface normal *along* the branch — leaves on horizontal branches went
/// edge-on to any horizontal viewer and effectively disappeared.
inline Quat quatLookDir(Vec3 forward) {
    Vec3 f = vnormOr(forward, {0, 1, 0});
    Vec3 up{0, 1, 0};
    // Surface normal: bisector of branch direction and world up. Falls back
    // to up when branch is straight up (degenerate bisector).
    Vec3 n = vnormOr(f + up, up);
    // First, rotate local +Y to the desired normal.
    Quat qN = quatFromTo(up, n);
    // Now spin around the new normal so local +Z aligns with branch tangent.
    Vec3 zAfter = quatRotate(qN, {0, 0, 1});
    Vec3 stemTarget = vnormOr(f - n * vdot(f, n), zAfter);
    Quat qS = quatFromTo(zAfter, stemTarget);
    return quatMul(qS, qN);
}

} // namespace plant_internal
} // namespace bromesh
