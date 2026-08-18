#include "bromesh/procedural/leaf_scatter.h"

#include "bromesh/manipulation/merge.h"
#include "bromesh/procedural/obstacle_field.h"

#include <bromath/bromath.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>

namespace bromesh {

using namespace bromath;

namespace {

constexpr float kTwoPi = bromath::TWO_PI;

// Fast per-segment RNG. Replaces std::mt19937_64 (2.5 KB of state, a heavy
// per-segment construction) in the parallel scatter loop: seeding is a single
// assignment and each draw is a handful of ops, so a segment that ends up
// placing no leaves costs almost nothing. splitmix64's statistical quality is
// ample for leaf jitter, and being purely functional in its state it stays
// deterministic per segment (same seed → same stream) under OpenMP.
struct FastRng {
    uint64_t s;
    explicit FastRng(uint64_t seed) : s(seed) {}
    inline uint64_t next() {
        uint64_t z = (s += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
    // 24-bit float in [0, 1).
    inline float uni01() { return (next() >> 40) * (1.0f / 16777216.0f); }
};

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
    // InstancedMeshNode canonical layout (16 floats per instance):
    // o[ 0..3 ] = row 0 of 3x4 affine matrix: r00, r01, r02, px
    // o[ 4..7 ] = row 1 of 3x4 affine matrix: r10, r11, r12, py
    // o[ 8..11] = row 2 of 3x4 affine matrix: r20, r21, r22, pz
    // o[12..15] = instance tint color (r, g, b, a) = (1, 1, 1, 1)
    out.push_back(sideAxis.x * scale);
    out.push_back(normal.x * scale);
    out.push_back(forward.x * scale);
    out.push_back(origin.x);

    out.push_back(sideAxis.y * scale);
    out.push_back(normal.y * scale);
    out.push_back(forward.y * scale);
    out.push_back(origin.y);

    out.push_back(sideAxis.z * scale);
    out.push_back(normal.z * scale);
    out.push_back(forward.z * scale);
    out.push_back(origin.z);

    out.push_back(1.0f);
    out.push_back(1.0f);
    out.push_back(1.0f);
    out.push_back(1.0f);
}

// Transform and stamp a prototype mesh at every placement in `pl` into a single merged MeshData.
MeshData stampMeshAtPlacements(const MeshData& srcMesh, const LeafPlacements& pl) {
    if (srcMesh.empty()) return {};
    const size_t count = pl.count();
    if (count == 0) return {};

    const size_t srcVCount = srcMesh.vertexCount();
    const size_t srcICount = srcMesh.indices.size();
    const bool hasNormals  = srcMesh.hasNormals();
    const bool hasUVs      = srcMesh.hasUVs();
    const bool hasColors   = srcMesh.hasColors();

    MeshData out;
    out.positions.resize(count * srcVCount * 3);
    if (hasNormals) out.normals.resize(count * srcVCount * 3);
    if (hasUVs)     out.uvs.resize(count * srcVCount * 2);
    if (hasColors)  out.colors.resize(count * srcVCount * 4);
    out.indices.resize(count * srcICount);

    for (size_t i = 0; i < count; ++i) {
        const float* M = &pl.transforms[i * 16];
        size_t vBase = i * srcVCount;
        uint32_t indexOffset = static_cast<uint32_t>(vBase);

        for (size_t v = 0; v < srcVCount; ++v) {
            size_t destPos = (vBase + v) * 3;
            float x = srcMesh.positions[v * 3 + 0];
            float y = srcMesh.positions[v * 3 + 1];
            float z = srcMesh.positions[v * 3 + 2];
            out.positions[destPos + 0] = M[0]*x + M[1]*y + M[2]*z  + M[3];
            out.positions[destPos + 1] = M[4]*x + M[5]*y + M[6]*z  + M[7];
            out.positions[destPos + 2] = M[8]*x + M[9]*y + M[10]*z + M[11];

            if (hasNormals) {
                float nx = srcMesh.normals[v * 3 + 0];
                float ny = srcMesh.normals[v * 3 + 1];
                float nz = srcMesh.normals[v * 3 + 2];
                float tx = M[0]*nx + M[1]*ny + M[2]*nz;
                float ty = M[4]*nx + M[5]*ny + M[6]*nz;
                float tz = M[8]*nx + M[9]*ny + M[10]*nz;
                float L = std::sqrt(tx*tx + ty*ty + tz*tz);
                if (L > 1e-8f) {
                    tx /= L; ty /= L; tz /= L;
                } else {
                    float nLen = std::sqrt(M[1]*M[1] + M[5]*M[5] + M[9]*M[9]);
                    if (nLen > 1e-8f) { tx = M[1]/nLen; ty = M[5]/nLen; tz = M[9]/nLen; }
                    else { tx = 0.0f; ty = 1.0f; tz = 0.0f; }
                }
                out.normals[destPos + 0] = tx;
                out.normals[destPos + 1] = ty;
                out.normals[destPos + 2] = tz;
            }

            if (hasUVs) {
                out.uvs[(vBase + v) * 2 + 0] = srcMesh.uvs[v * 2 + 0];
                out.uvs[(vBase + v) * 2 + 1] = srcMesh.uvs[v * 2 + 1];
            }
            if (hasColors) {
                out.colors[(vBase + v) * 4 + 0] = srcMesh.colors[v * 4 + 0];
                out.colors[(vBase + v) * 4 + 1] = srcMesh.colors[v * 4 + 1];
                out.colors[(vBase + v) * 4 + 2] = srcMesh.colors[v * 4 + 2];
                out.colors[(vBase + v) * 4 + 3] = srcMesh.colors[v * 4 + 3];
            }
        }

        size_t iBase = i * srcICount;
        for (size_t idx = 0; idx < srcICount; ++idx) {
            out.indices[iBase + idx] = srcMesh.indices[idx] + indexOffset;
        }
    }

    return out;
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

    const Vec3 worldUp{0, 1, 0};
    const float refRadius = (opts.maxRadius > 1e-8f) ? opts.maxRadius : 1.0f;
    const size_t segCount = segments.size();

    std::vector<std::vector<float>> segTransforms(segCount);
    std::vector<std::vector<float>> segRadius(segCount);
    std::vector<std::vector<int>> segDepth(segCount);

    #pragma omp parallel for schedule(dynamic, 16) if(segCount > 32)
    for (int i = 0; i < static_cast<int>(segCount); ++i) {
        const BranchSegment& seg = segments[static_cast<size_t>(i)];
        Vec3 d = seg.to - seg.from;
        float length = vlen(d);
        if (length < 1e-6f) continue;
        if (seg.depth < opts.minDepth) continue;
        if (seg.radius > 0.0f && seg.radius > opts.maxRadius) continue;
        if (opts.terminalOnly && childCount[static_cast<size_t>(i)] > 0) continue;

        Vec3 T = d * (1.0f / length);

        float weight = 1.0f;
        if (!opts.densityWeight.empty()) {
            weight = (static_cast<size_t>(i) < opts.densityWeight.size())
                         ? std::max(0.0f, opts.densityWeight[static_cast<size_t>(i)])
                         : 0.0f;
            if (weight <= 0.0f) continue;
        }

        FastRng rng(opts.seed ^ (static_cast<uint64_t>(i) * 0x9E3779B97F4A7C15ULL));

        float expected = length * opts.perUnitLength * weight;
        int sampleCount = static_cast<int>(std::floor(expected));
        if (rng.uni01() < (expected - static_cast<float>(sampleCount))) {
            ++sampleCount;
        }
        if (sampleCount <= 0) continue;

        auto& transforms = segTransforms[static_cast<size_t>(i)];
        auto& radii = segRadius[static_cast<size_t>(i)];
        auto& depths = segDepth[static_cast<size_t>(i)];
        transforms.reserve(static_cast<size_t>(sampleCount) * 16);
        radii.reserve(static_cast<size_t>(sampleCount));
        depths.reserve(static_cast<size_t>(sampleCount));

        for (int s = 0; s < sampleCount; ++s) {
            float u = rng.uni01();
            float t = (opts.densityFalloff > 0.0f)
                ? 1.0f - std::pow(1.0f - u, 1.0f + opts.densityFalloff)
                : u;
            Vec3 P = seg.from + d * t;

            // Obstacle test on the candidate origin.
            const int selfTag = static_cast<int>(i);
            if (opts.avoid != nullptr && !opts.avoid->empty()) {
                if (opts.avoid->tooClose(P, opts.obstacleClearance, selfTag)) {
                    if (opts.obstaclePushout > 0.0f) {
                        auto n = opts.avoid->nearest(P, selfTag);
                        if (n.tag != -1 || std::isfinite(n.distance)) {
                            P = n.point + n.normal * (opts.obstacleClearance + opts.obstaclePushout);
                        }
                        if (opts.avoid->tooClose(P, opts.obstacleClearance, selfTag)) {
                            continue;
                        }
                    } else {
                        continue;
                    }
                }
            }
            if (!opts.keepOut.empty()) {
                bool blocked = false;
                for (const Sphere& sp : opts.keepOut) {
                    Vec3 dv = P - sp.center;
                    float r = sp.radius + opts.obstacleClearance;
                    if (vdot(dv, dv) <= r * r) { blocked = true; break; }
                }
                if (blocked) continue;
            }

            float phi = rng.uni01() * kTwoPi;
            Vec3 e1 = perpendicularUnit(T);
            Vec3 e2 = vcross(T, e1);
            Vec3 R = e1 * std::cos(phi) + e2 * std::sin(phi);

            Vec3 Fraw = R * (1.0f - opts.upBias) + worldUp * opts.upBias;
            Vec3 F = vnormOr(Fraw, R);

            if (opts.tiltJitter > 0.0f) {
                float tiltAng = (rng.uni01() * 2.0f - 1.0f) * opts.tiltJitter;
                Vec3 tiltAxis = vnormOr(vcross(F, T), e2);
                F = vnorm(qrotate(qaxisAngle(tiltAxis, tiltAng), F));
            }

            Vec3 sideAxis = vcross(F, worldUp);
            if (vdot(sideAxis, sideAxis) < 1e-8f) {
                sideAxis = vcross(F, Vec3{1, 0, 0});
            }
            sideAxis = vnorm(sideAxis);
            Vec3 N = vnorm(vcross(sideAxis, F));

            if (opts.rollJitter > 0.0f) {
                float rollAng = (rng.uni01() * 2.0f - 1.0f) * opts.rollJitter;
                Quat q = qaxisAngle(F, rollAng);
                sideAxis = vnorm(qrotate(q, sideAxis));
                N = vnorm(qrotate(q, N));
            }

            float jitter = (rng.uni01() * 2.0f - 1.0f) * opts.scaleJitter;
            float radiusFactor = 1.0f;
            if (opts.scaleByRadius > 0.0f) {
                float ratio = (seg.radius > 0.0f) ? (seg.radius / refRadius) : 1.0f;
                radiusFactor = 1.0f + opts.scaleByRadius * (ratio - 1.0f);
                if (radiusFactor < 0.05f) radiusFactor = 0.05f;
            }
            float scale = opts.baseScale * (1.0f + jitter) * radiusFactor;
            if (scale < 1e-6f) continue;

            writeMatrix(transforms, sideAxis, N, F, P, scale);
            radii.push_back(seg.radius);
            depths.push_back(seg.depth);
        }
    }

    size_t totalTransforms = 0;
    for (size_t i = 0; i < segCount; ++i) {
        totalTransforms += segRadius[i].size();
    }
    out.transforms.reserve(totalTransforms * 16);
    out.branchRadius.reserve(totalTransforms);
    out.branchDepth.reserve(totalTransforms);

    if (opts.dedupRadius > 0.0f) {
        SpatialHash3D dedupHash(opts.dedupRadius);
        int32_t dedupId = 0;
        std::vector<int32_t> nearby;
        for (size_t i = 0; i < segCount; ++i) {
            const auto& transforms = segTransforms[i];
            const size_t placed = segRadius[i].size();
            for (size_t k = 0; k < placed; ++k) {
                const float* M = transforms.data() + k * 16;
                Vec3 P{M[3], M[7], M[11]};
                nearby.clear();
                dedupHash.radiusQuery(P, opts.dedupRadius, nearby);
                if (!nearby.empty()) continue;
                dedupHash.insert(P, dedupId++);
                out.transforms.insert(out.transforms.end(), M, M + 16);
                out.branchRadius.push_back(segRadius[i][k]);
                out.branchDepth.push_back(segDepth[i][k]);
            }
        }
    } else {
        for (size_t i = 0; i < segCount; ++i) {
            out.transforms.insert(out.transforms.end(), segTransforms[i].begin(), segTransforms[i].end());
            out.branchRadius.insert(out.branchRadius.end(), segRadius[i].begin(), segRadius[i].end());
            out.branchDepth.insert(out.branchDepth.end(), segDepth[i].begin(), segDepth[i].end());
        }
    }

    return out;
}

LeafPlacements placeLeafClustersOnBranches(
    const std::vector<BranchSegment>& segments,
    const LeafPlacementOptions& opts) {
    LeafPlacements out;
    if (segments.empty()) return out;

    std::vector<int> childCount(segments.size(), 0);
    for (size_t i = 0; i < segments.size(); ++i) {
        int p = segments[i].parent;
        if (p >= 0 && static_cast<size_t>(p) < segments.size()) {
            ++childCount[p];
        }
    }

    const Vec3 worldUp{0, 1, 0};
    const Vec3 gravity{0, -1, 0};
    const float refRadius = (opts.maxRadius > 1e-8f) ? opts.maxRadius : 1.0f;
    const size_t segCount = segments.size();

    std::vector<std::vector<float>> segTransforms(segCount);
    std::vector<std::vector<float>> segRadius(segCount);
    std::vector<std::vector<int>> segDepth(segCount);

    #pragma omp parallel for schedule(dynamic, 16) if(segCount > 32)
    for (int i = 0; i < static_cast<int>(segCount); ++i) {
        const BranchSegment& seg = segments[static_cast<size_t>(i)];
        Vec3 d = seg.to - seg.from;
        float length = vlen(d);
        if (length < 1e-6f) continue;
        if (seg.depth < opts.minDepth) continue;
        if (seg.radius > 0.0f && seg.radius > opts.maxRadius) continue;
        if (opts.terminalOnly && childCount[static_cast<size_t>(i)] > 0) continue;

        Vec3 T = d * (1.0f / length);

        float weight = 1.0f;
        if (!opts.densityWeight.empty()) {
            weight = (static_cast<size_t>(i) < opts.densityWeight.size())
                         ? std::max(0.0f, opts.densityWeight[static_cast<size_t>(i)])
                         : 0.0f;
            if (weight <= 0.0f) continue;
        }

        FastRng rng(opts.seed ^ (static_cast<uint64_t>(i) * 0x9E3779B97F4A7C15ULL));

        float expected = length * opts.perUnitLength * weight;
        int sampleCount = static_cast<int>(std::floor(expected));
        if (rng.uni01() < (expected - static_cast<float>(sampleCount))) {
            ++sampleCount;
        }
        if (sampleCount <= 0) continue;

        auto& transforms = segTransforms[static_cast<size_t>(i)];
        auto& radii = segRadius[static_cast<size_t>(i)];
        auto& depths = segDepth[static_cast<size_t>(i)];
        transforms.reserve(static_cast<size_t>(sampleCount) * 16);
        radii.reserve(static_cast<size_t>(sampleCount));
        depths.reserve(static_cast<size_t>(sampleCount));

        for (int s = 0; s < sampleCount; ++s) {
            float u = rng.uni01();
            float t = (opts.densityFalloff > 0.0f)
                ? 1.0f - std::pow(1.0f - u, 1.0f + opts.densityFalloff)
                : u;

            float phi = rng.uni01() * kTwoPi;
            Vec3 e1 = perpendicularUnit(T);
            Vec3 e2 = vcross(T, e1);
            Vec3 R = e1 * std::cos(phi) + e2 * std::sin(phi);

            // Attachment point on branch surface
            Vec3 P = seg.from + d * t + R * seg.radius;

            // Obstacle test on candidate origin
            const int selfTag = static_cast<int>(i);
            if (opts.avoid != nullptr && !opts.avoid->empty()) {
                if (opts.avoid->tooClose(P, opts.obstacleClearance, selfTag)) {
                    if (opts.obstaclePushout > 0.0f) {
                        auto n = opts.avoid->nearest(P, selfTag);
                        if (n.tag != -1 || std::isfinite(n.distance)) {
                            P = n.point + n.normal * (opts.obstacleClearance + opts.obstaclePushout);
                        }
                        if (opts.avoid->tooClose(P, opts.obstacleClearance, selfTag)) {
                            continue;
                        }
                    } else {
                        continue;
                    }
                }
            }
            if (!opts.keepOut.empty()) {
                bool blocked = false;
                for (const Sphere& sp : opts.keepOut) {
                    Vec3 dv = P - sp.center;
                    float r = sp.radius + opts.obstacleClearance;
                    if (vdot(dv, dv) <= r * r) { blocked = true; break; }
                }
                if (blocked) continue;
            }

            // Cluster growth direction: along twig growth direction (+tangent),
            // flared slightly outward laterally on non-terminal points
            float lateralFlare = (t > 0.85f && childCount[static_cast<size_t>(i)] == 0) ? 0.15f : 0.45f;
            Vec3 D = vnormOr(T * (1.0f - lateralFlare) + R * lateralFlare, T);

            // Gravity sag / droop
            float horiz = std::sqrt(D.x * D.x + D.z * D.z);
            Vec3 F = vnormOr(D + gravity * (0.15f * horiz), D);

            // Phototropism: turn cluster toward light
            Vec3 Fraw = F * (1.0f - opts.upBias * 0.4f) + worldUp * (opts.upBias * 0.4f);
            F = vnormOr(Fraw, F);

            // Tilt jitter
            if (opts.tiltJitter > 0.0f) {
                float tiltAng = (rng.uni01() * 2.0f - 1.0f) * opts.tiltJitter;
                Vec3 tiltAxis = vnormOr(vcross(F, T), e2);
                F = vnorm(qrotate(qaxisAngle(tiltAxis, tiltAng), F));
            }

            // Upper surface normal faces skyward / light
            Vec3 sideAxis = vcross(F, worldUp);
            if (vdot(sideAxis, sideAxis) < 1e-8f) {
                sideAxis = vcross(F, Vec3{1, 0, 0});
            }
            sideAxis = vnorm(sideAxis);
            Vec3 N = vnorm(vcross(sideAxis, F));

            // Roll jitter
            if (opts.rollJitter > 0.0f) {
                float rollAng = (rng.uni01() * 2.0f - 1.0f) * opts.rollJitter;
                Quat q = qaxisAngle(F, rollAng);
                sideAxis = vnorm(qrotate(q, sideAxis));
                N = vnorm(qrotate(q, N));
            }

            float jitter = (rng.uni01() * 2.0f - 1.0f) * opts.scaleJitter;
            float radiusFactor = 1.0f;
            if (opts.scaleByRadius > 0.0f) {
                float ratio = (seg.radius > 0.0f) ? (seg.radius / refRadius) : 1.0f;
                radiusFactor = 1.0f + opts.scaleByRadius * (ratio - 1.0f);
                if (radiusFactor < 0.05f) radiusFactor = 0.05f;
            }
            float scale = opts.baseScale * (1.0f + jitter) * radiusFactor;
            if (scale < 1e-6f) continue;

            writeMatrix(transforms, sideAxis, N, F, P, scale);
            radii.push_back(seg.radius);
            depths.push_back(seg.depth);
        }
    }

    size_t totalTransforms = 0;
    for (size_t i = 0; i < segCount; ++i) {
        totalTransforms += segRadius[i].size();
    }
    out.transforms.reserve(totalTransforms * 16);
    out.branchRadius.reserve(totalTransforms);
    out.branchDepth.reserve(totalTransforms);

    if (opts.dedupRadius > 0.0f) {
        SpatialHash3D dedupHash(opts.dedupRadius);
        int32_t dedupId = 0;
        std::vector<int32_t> nearby;
        for (size_t i = 0; i < segCount; ++i) {
            const auto& transforms = segTransforms[i];
            const size_t placed = segRadius[i].size();
            for (size_t k = 0; k < placed; ++k) {
                const float* M = transforms.data() + k * 16;
                Vec3 P{M[3], M[7], M[11]};
                nearby.clear();
                dedupHash.radiusQuery(P, opts.dedupRadius, nearby);
                if (!nearby.empty()) continue;
                dedupHash.insert(P, dedupId++);
                out.transforms.insert(out.transforms.end(), M, M + 16);
                out.branchRadius.push_back(segRadius[i][k]);
                out.branchDepth.push_back(segDepth[i][k]);
            }
        }
    } else {
        for (size_t i = 0; i < segCount; ++i) {
            out.transforms.insert(out.transforms.end(), segTransforms[i].begin(), segTransforms[i].end());
            out.branchRadius.insert(out.branchRadius.end(), segRadius[i].begin(), segRadius[i].end());
            out.branchDepth.insert(out.branchDepth.end(), segDepth[i].begin(), segDepth[i].end());
        }
    }

    return out;
}

MeshData scatterLeaves(
    const std::vector<BranchSegment>& segments,
    const MeshData& leaf,
    const LeafPlacementOptions& opts) {
    if (leaf.empty()) return {};
    LeafPlacements pl = placeLeavesOnBranches(segments, opts);
    return stampMeshAtPlacements(leaf, pl);
}

MeshData scatterLeafClusters(
    const std::vector<BranchSegment>& segments,
    Phyllotaxy phyllotaxy,
    const LeafClusterOptions& clusterOpts,
    const LeafPlacementOptions& scatterOpts) {
    MeshData cluster = leafCluster(phyllotaxy, clusterOpts);
    if (cluster.empty()) return {};
    LeafPlacements pl = placeLeafClustersOnBranches(segments, scatterOpts);
    return stampMeshAtPlacements(cluster, pl);
}

std::vector<int> packAnchors(
    const std::vector<Vec3>& candidates,
    const CapsuleField* avoid,
    const std::vector<Sphere>& reservedKeepOut,
    const AnchorPackOptions& opts) {
    std::vector<int> accepted;
    if (candidates.empty()) return accepted;

    std::vector<int> order(candidates.size());
    std::iota(order.begin(), order.end(), 0);
    std::mt19937_64 rng(opts.seed);
    std::shuffle(order.begin(), order.end(), rng);

    const bool useSpacing = opts.minSpacing > 0.0f;
    SpatialHash3D spacingHash(useSpacing ? opts.minSpacing : 1.0f);
    std::vector<int32_t> nearby;

    const int cap = (opts.maxCount > 0) ? opts.maxCount : -1;

    for (int idx : order) {
        if (cap >= 0 && static_cast<int>(accepted.size()) >= cap) break;
        Vec3 p = candidates[static_cast<size_t>(idx)];

        if (avoid != nullptr && !avoid->empty()) {
            if (avoid->tooClose(p, opts.minObstacleDistance, -1)) continue;
        }
        if (!reservedKeepOut.empty()) {
            bool blocked = false;
            for (const Sphere& sp : reservedKeepOut) {
                Vec3 dv = p - sp.center;
                float r = sp.radius;
                if (vdot(dv, dv) <= r * r) { blocked = true; break; }
            }
            if (blocked) continue;
        }
        if (useSpacing) {
            nearby.clear();
            spacingHash.radiusQuery(p, opts.minSpacing, nearby);
            if (!nearby.empty()) continue;
            spacingHash.insert(p, static_cast<int32_t>(accepted.size()));
        }

        accepted.push_back(idx);
    }
    return accepted;
}

} // namespace bromesh
