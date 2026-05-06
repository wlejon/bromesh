#pragma once

#include "bromesh/manipulation/sweep.h"
#include "bromesh/mesh_data.h"
#include "bromesh/procedural/vec_math.h"

#include <vector>

namespace bromesh {

struct BezierSweepOptions {
    /// Total number of path samples along the spline (>= 4).
    int samples = 32;
    /// Per-sample profile scale (size 0 or 1 = constant; otherwise
    /// length must be `samples`, evaluated at t = i / (samples - 1)).
    std::vector<float> profileScale;
    /// Per-sample twist in radians (same shape rules as profileScale).
    std::vector<float> twist;
    bool capStart = true;
    bool capEnd = true;
    bool closeProfile = true;
    bool miterJoints = true;
};

/// Sweep a 2D `profile` along a cubic-bezier polyline. `controlPoints` is
/// interpreted as a sequence of cubic bezier segments where every group of
/// 4 points (P0,P1,P2,P3) defines one segment, and consecutive segments
/// SHARE their endpoint — i.e. for two segments you pass 7 points
/// (P0,P1,P2,P3=Q0,Q1,Q2,Q3). The total length of `controlPoints` must
/// satisfy `(N - 1) % 3 == 0` and `N >= 4`. Empty / malformed input
/// returns an empty MeshData.
MeshData bezierSweep(const std::vector<Vec3>& controlPoints,
                     const std::vector<Vec2>& profile,
                     const BezierSweepOptions& opts = {});

} // namespace bromesh
