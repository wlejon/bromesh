#include "bromesh/procedural/leaf_clumps.h"

#include <bromath/bromath.h>

#include <algorithm>
#include <cmath>

namespace bromesh {

using namespace bromath;

namespace {

constexpr float kTwoPi = bromath::TWO_PI;
constexpr uint32_t kGolden = 0x9e3779b9u;

// Draw `k` of a stream keyed on `base`, in [0, 1).
inline float draw(uint32_t base, uint32_t k) {
    return static_cast<float>(leafClumpHash(base + k * kGolden) >> 8) * (1.0f / 16777216.0f);
}

inline uint32_t clumpBase(uint32_t index, uint32_t seed) {
    return leafClumpHash(seed ^ leafClumpHash(index * 0x632be5abu + 0x51ed27u));
}

Vec3 perpendicularUnit(Vec3 t) {
    Vec3 c = vcross(t, Vec3{0.0f, 1.0f, 0.0f});
    if (vdot(c, c) < 1e-8f) c = vcross(t, Vec3{1.0f, 0.0f, 0.0f});
    return vnorm(c);
}

}  // namespace

void leafClumpLobe(const LeafClump& c, uint32_t index, int lobe, uint32_t seed,
                   Vec3& center, float& radius) {
    if (lobe <= 0) {
        center = c.center;
        radius = c.radius * 0.72f;
        return;
    }
    const uint32_t b = leafClumpHash(clumpBase(index, seed) + static_cast<uint32_t>(lobe) * kGolden);
    // Lobes spread sideways and up more than down: clumps sit like cumulus.
    const float z = draw(b, 0) * 1.2f - 0.4f;
    const float phi = draw(b, 1) * kTwoPi;
    const float s = std::sqrt(std::max(0.0f, 1.0f - z * z));
    const Vec3 d{s * std::cos(phi), z, s * std::sin(phi)};
    center = c.center + d * (c.radius * 0.42f);
    radius = c.radius * (0.5f + 0.2f * draw(b, 2));
}

int leafClumpCount(const LeafClump& c, uint32_t index, const LeafClumpOptions& opts) {
    const float area = 4.0f * bromath::PI * c.radius * c.radius;
    float expected = opts.perArea * std::max(0.0f, c.density) * area;
    if (!(expected > 0.0f)) return 0;
    expected = std::min(expected, static_cast<float>(kMaxLeavesPerSegment));
    const float frac = draw(clumpBase(index, opts.seed), 7);
    return static_cast<int>(std::floor(expected + frac));
}

bool leafInClump(const LeafClump& c, uint32_t index, uint32_t leaf, const LeafClumpOptions& opts,
                 float o[16]) {
    const uint32_t base = leafClumpHash(opts.seed ^ leafClumpHash(leaf + 1u));
    const int lobe = std::min(kLeafClumpLobes - 1, static_cast<int>(draw(base, 0) * kLeafClumpLobes));
    Vec3 lc;
    float lr;
    leafClumpLobe(c, index, lobe, opts.seed, lc, lr);

    // A point in the lobe's outer shell, the clump flattened underneath.
    const float z = draw(base, 1) * 2.0f - 1.0f;
    const float phi = draw(base, 2) * kTwoPi;
    const float s = std::sqrt(std::max(0.0f, 1.0f - z * z));
    const Vec3 u{s * std::cos(phi), z, s * std::sin(phi)};
    Vec3 off = u * (lr * (1.0f - opts.shell * draw(base, 3)));
    off.y *= u.y < 0.0f ? 0.7f : 0.9f;
    const Vec3 P = lc + off;

    // Facing out of the lobe, leaning up, jittered.
    Vec3 N = vnorm(u * (1.0f - opts.upBias) + Vec3{0.0f, opts.upBias, 0.0f});
    const Vec3 jit{draw(base, 4) * 2.0f - 1.0f, draw(base, 5) * 2.0f - 1.0f, draw(base, 6) * 2.0f - 1.0f};
    N = N + jit * opts.tiltJitter;
    if (vdot(N, N) < 1e-8f) N = u;
    N = vnorm(N);

    // The tip any way round the normal.
    const float theta = draw(base, 7) * kTwoPi;
    const Vec3 t0 = perpendicularUnit(N);
    const Vec3 t1 = vcross(N, t0);
    const Vec3 F = vnorm(t0 * std::cos(theta) + t1 * std::sin(theta));
    const Vec3 X = vcross(F, N);

    const float scale = opts.baseScale * (1.0f + (draw(base, 8) * 2.0f - 1.0f) * opts.scaleJitter);
    if (!(scale > 1e-6f)) return false;
    o[0] = X.x * scale; o[1] = N.x * scale; o[2] = F.x * scale;  o[3] = P.x;
    o[4] = X.y * scale; o[5] = N.y * scale; o[6] = F.y * scale;  o[7] = P.y;
    o[8] = X.z * scale; o[9] = N.z * scale; o[10] = F.z * scale; o[11] = P.z;
    o[12] = o[13] = o[14] = o[15] = 1.0f;
    return true;
}

LeafPlacements placeLeavesInClumps(const std::vector<LeafClump>& clumps, const LeafClumpOptions& opts) {
    LeafPlacements out;
    uint32_t leaf = 0;
    float rec[16];
    for (size_t ci = 0; ci < clumps.size(); ++ci) {
        const int n = leafClumpCount(clumps[ci], static_cast<uint32_t>(ci), opts);
        for (int k = 0; k < n; ++k, ++leaf) {
            if (!leafInClump(clumps[ci], static_cast<uint32_t>(ci), leaf, opts, rec)) continue;
            out.transforms.insert(out.transforms.end(), rec, rec + 16);
            out.branchRadius.push_back(clumps[ci].radius);
            out.branchDepth.push_back(0);
            out.segment.push_back(static_cast<int32_t>(ci));
        }
    }
    return out;
}

} // namespace bromesh
