#include "bromesh/manipulation/bezier_sweep.h"

#include <bromath/bromath.h>

#include <algorithm>
#include <cmath>

namespace bromesh {

using namespace bromath;

namespace {

void resample(const std::vector<float>& src, std::vector<float>& dst, int n) {
    dst.clear();
    if (src.empty() || n <= 0) return;
    if (src.size() == 1) { dst.assign(n, src[0]); return; }
    if (static_cast<int>(src.size()) == n) { dst = src; return; }
    dst.reserve(n);
    int last = static_cast<int>(src.size()) - 1;
    for (int i = 0; i < n; ++i) {
        float t = (n == 1) ? 0.0f : static_cast<float>(i) / static_cast<float>(n - 1);
        float pos = t * static_cast<float>(last);
        int i0 = static_cast<int>(std::floor(pos));
        int i1 = std::min(last, i0 + 1);
        float f = pos - static_cast<float>(i0);
        dst.push_back(src[i0] * (1.0f - f) + src[i1] * f);
    }
}

} // namespace

MeshData bezierSweep(const std::vector<Vec3>& controlPoints,
                     const std::vector<Vec2>& profile,
                     const BezierSweepOptions& opts) {
    int N = static_cast<int>(controlPoints.size());
    if (N < 4 || ((N - 1) % 3) != 0) return {};
    int segCount = (N - 1) / 3;
    int samples = std::max(4, opts.samples);

    std::vector<Vec3> path;
    path.reserve(samples);
    // Sample globally over t in [0,1], distributed evenly across the multi-segment spline.
    for (int i = 0; i < samples; ++i) {
        float t = static_cast<float>(i) / static_cast<float>(samples - 1);
        // Map t to a segment + local s.
        float scaled = t * static_cast<float>(segCount);
        int seg = std::min(segCount - 1, static_cast<int>(std::floor(scaled)));
        float s = scaled - static_cast<float>(seg);
        int base = seg * 3;
        Vec3 pt = cbezier(
            controlPoints[base + 0],
            controlPoints[base + 1],
            controlPoints[base + 2],
            controlPoints[base + 3],
            s);
        path.push_back(pt);
    }

    SweepOptions sw;
    sw.capStart = opts.capStart;
    sw.capEnd = opts.capEnd;
    sw.closeProfile = opts.closeProfile;
    sw.miterJoints = opts.miterJoints;
    resample(opts.profileScale, sw.profileScale, samples);
    resample(opts.twist, sw.twist, samples);

    return sweep(profile, path, sw);
}

} // namespace bromesh
