#include "bromesh/procedural/leaf_scatter.h"

#include "bromesh/manipulation/merge.h"
#include "bromesh/optimization/spatial_hash.h"
#include "bromesh/procedural/vec_math.h"

#include <algorithm>
#include <cmath>
#include <random>

namespace bromesh {

namespace {

constexpr float kPi = 3.14159265358979323846f;
constexpr float kTwoPi = 2.0f * kPi;

// Pick any unit vector perpendicular to `t` (assumed unit length).
Vec3 perpendicularUnit(Vec3 t) {
    Vec3 worldUp{0, 1, 0};
    Vec3 c = vcross(t, worldUp);
    if (vdot(c, c) < 1e-8f) {
        c = vcross(t, Vec3{1, 0, 0});
    }
    return vnorm(c);
}

void writeMatrix(std::vector<float>& out,
                 Vec3 sideAxis, Vec3 normal, Vec3 forward,
                 Vec3 origin, float scale) {
    // Column-major 4x4: cols are X-basis, Y-basis, Z-basis, translation.
    out.push_back(sideAxis.x * scale);
    out.push_back(sideAxis.y * scale);
    out.push_back(sideAxis.z * scale);
    out.push_back(0.0f);

    out.push_back(normal.x * scale);
    out.push_back(normal.y * scale);
    out.push_back(normal.z * scale);
    out.push_back(0.0f);

    out.push_back(forward.x * scale);
    out.push_back(forward.y * scale);
    out.push_back(forward.z * scale);
    out.push_back(0.0f);

    out.push_back(origin.x);
    out.push_back(origin.y);
    out.push_back(origin.z);
    out.push_back(1.0f);
}

} // namespace

LeafPlacements placeLeavesOnBranches(
    const std::vector<BranchSegment>& segments,
    const LeafPlacementOptions& opts) {
    LeafPlacements out;
    if (segments.empty()) return out;

    // Child counts for terminalOnly filter.
    std::vector<int> childCount(segments.size(), 0);
    for (size_t i = 0; i < segments.size(); ++i) {
        int p = segments[i].parent;
        if (p >= 0 && static_cast<size_t>(p) < segments.size()) {
            ++childCount[p];
        }
    }

    std::mt19937_64 rng(opts.seed);
    std::uniform_real_distribution<float> uni01(0.0f, 1.0f);

    SpatialHash3D dedupHash(opts.dedupRadius > 0.0f
                                ? opts.dedupRadius
                                : 1.0f);
    int32_t dedupId = 0;
    std::vector<int32_t> nearby;

    const Vec3 worldUp{0, 1, 0};
    const float refRadius = (opts.maxRadius > 1e-8f) ? opts.maxRadius : 1.0f;

    for (size_t i = 0; i < segments.size(); ++i) {
        const BranchSegment& seg = segments[i];
        Vec3 d = seg.to - seg.from;
        float length = vlen(d);
        if (length < 1e-6f) continue;
        if (seg.depth < opts.minDepth) continue;
        if (seg.radius > 0.0f && seg.radius > opts.maxRadius) continue;
        if (opts.terminalOnly && childCount[i] > 0) continue;

        Vec3 T = d * (1.0f / length);

        float expected = length * opts.perUnitLength;
        // Stochastic rounding so very short segments still occasionally place.
        int sampleCount = static_cast<int>(std::floor(expected));
        if (uni01(rng) < (expected - static_cast<float>(sampleCount))) {
            ++sampleCount;
        }
        if (sampleCount <= 0) continue;

        for (int s = 0; s < sampleCount; ++s) {
            float u = uni01(rng);
            float t = (opts.densityFalloff > 0.0f)
                ? 1.0f - std::pow(1.0f - u, 1.0f + opts.densityFalloff)
                : u;
            Vec3 P = seg.from + d * t;

            // Azimuth around the branch.
            float phi = uni01(rng) * kTwoPi;
            Vec3 e1 = perpendicularUnit(T);
            Vec3 e2 = vcross(T, e1);
            Vec3 R = e1 * std::cos(phi) + e2 * std::sin(phi);

            // Forward = blend(radial, world up).
            Vec3 Fraw = R * (1.0f - opts.upBias) + worldUp * opts.upBias;
            Vec3 F = vnormOr(Fraw, R);

            // Tilt: pitch F forward/back along T using axis perpendicular to F and T.
            if (opts.tiltJitter > 0.0f) {
                float tiltAng = (uni01(rng) * 2.0f - 1.0f) * opts.tiltJitter;
                Vec3 tiltAxis = vnormOr(vcross(F, T), e2);
                F = vnorm(quatRotate(quatAxisAngle(tiltAxis, tiltAng), F));
            }

            // Side / normal — N points roughly upward when F is roughly horizontal.
            Vec3 sideAxis = vcross(F, worldUp);
            if (vdot(sideAxis, sideAxis) < 1e-8f) {
                sideAxis = vcross(F, Vec3{1, 0, 0});
            }
            sideAxis = vnorm(sideAxis);
            Vec3 N = vnorm(vcross(sideAxis, F));

            // Roll around F.
            if (opts.rollJitter > 0.0f) {
                float rollAng = (uni01(rng) * 2.0f - 1.0f) * opts.rollJitter;
                Quat q = quatAxisAngle(F, rollAng);
                sideAxis = vnorm(quatRotate(q, sideAxis));
                N = vnorm(quatRotate(q, N));
            }

            // Scale.
            float jitter = (uni01(rng) * 2.0f - 1.0f) * opts.scaleJitter;
            float radiusFactor = 1.0f;
            if (opts.scaleByRadius > 0.0f) {
                float ratio = (seg.radius > 0.0f) ? (seg.radius / refRadius) : 1.0f;
                radiusFactor = 1.0f + opts.scaleByRadius * (ratio - 1.0f);
                if (radiusFactor < 0.05f) radiusFactor = 0.05f;
            }
            float scale = opts.baseScale * (1.0f + jitter) * radiusFactor;
            if (scale < 1e-6f) continue;

            // Dedup.
            if (opts.dedupRadius > 0.0f) {
                nearby.clear();
                dedupHash.radiusQuery(P, opts.dedupRadius, nearby);
                if (!nearby.empty()) continue;
                dedupHash.insert(P, dedupId++);
            }

            writeMatrix(out.transforms, sideAxis, N, F, P, scale);
            out.branchRadius.push_back(seg.radius);
            out.branchDepth.push_back(seg.depth);
        }
    }

    return out;
}

namespace {

void transformLeafCopy(const MeshData& src, const float* M, MeshData& dst) {
    size_t vcount = src.vertexCount();
    dst.positions.resize(vcount * 3);
    for (size_t i = 0; i < vcount; ++i) {
        float x = src.positions[i * 3 + 0];
        float y = src.positions[i * 3 + 1];
        float z = src.positions[i * 3 + 2];
        dst.positions[i * 3 + 0] = M[0]*x + M[4]*y + M[8]*z  + M[12];
        dst.positions[i * 3 + 1] = M[1]*x + M[5]*y + M[9]*z  + M[13];
        dst.positions[i * 3 + 2] = M[2]*x + M[6]*y + M[10]*z + M[14];
    }
    if (src.hasNormals()) {
        dst.normals.resize(vcount * 3);
        for (size_t i = 0; i < vcount; ++i) {
            float x = src.normals[i * 3 + 0];
            float y = src.normals[i * 3 + 1];
            float z = src.normals[i * 3 + 2];
            float nx = M[0]*x + M[4]*y + M[8]*z;
            float ny = M[1]*x + M[5]*y + M[9]*z;
            float nz = M[2]*x + M[6]*y + M[10]*z;
            float L = std::sqrt(nx*nx + ny*ny + nz*nz);
            if (L > 1e-8f) { nx /= L; ny /= L; nz /= L; }
            dst.normals[i * 3 + 0] = nx;
            dst.normals[i * 3 + 1] = ny;
            dst.normals[i * 3 + 2] = nz;
        }
    }
    dst.uvs = src.uvs;
    dst.colors = src.colors;
    dst.indices = src.indices;
}

} // namespace

MeshData scatterLeaves(
    const std::vector<BranchSegment>& segments,
    const MeshData& leaf,
    const LeafPlacementOptions& opts) {
    if (leaf.empty()) return {};
    LeafPlacements pl = placeLeavesOnBranches(segments, opts);
    if (pl.count() == 0) return {};

    std::vector<MeshData> parts;
    parts.reserve(pl.count());
    MeshData stamped;
    for (size_t i = 0; i < pl.count(); ++i) {
        transformLeafCopy(leaf, &pl.transforms[i * 16], stamped);
        parts.push_back(stamped);
    }
    return mergeMeshes(parts);
}

} // namespace bromesh
