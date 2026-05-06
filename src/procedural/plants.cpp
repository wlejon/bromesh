#include "bromesh/procedural/plants.h"

#include "bromesh/manipulation/merge.h"
#include "bromesh/manipulation/sweep.h"
#include "bromesh/manipulation/transform.h"
#include "bromesh/primitives/primitives.h"
#include "bromesh/procedural/vec_math.h"

#include <algorithm>
#include <cmath>

namespace bromesh {

namespace {

constexpr int kAtlasCols = 4;
constexpr int kAtlasRows = 4;

void atlasCell(LeafShape shape, float& u0, float& v0, float& u1, float& v1) {
    int idx = static_cast<int>(shape);
    int col = idx % kAtlasCols;
    int row = idx / kAtlasCols;
    float du = 1.0f / kAtlasCols;
    float dv = 1.0f / kAtlasRows;
    u0 = col * du;
    v0 = row * dv;
    u1 = u0 + du;
    v1 = v0 + dv;
}

} // namespace

MeshData leafCard(LeafShape shape, const LeafCardOptions& opts) {
    int wseg = std::max(1, opts.widthSegments);
    int lseg = std::max(1, opts.lengthSegments);
    int wcount = wseg + 1;
    int lcount = lseg + 1;

    float u0, v0, u1, v1;
    if (opts.fullUV) {
        u0 = 0.0f; v0 = 0.0f; u1 = 1.0f; v1 = 1.0f;
    } else {
        atlasCell(shape, u0, v0, u1, v1);
    }

    MeshData m;
    m.positions.reserve(wcount * lcount * 3);
    m.normals.reserve(wcount * lcount * 3);
    m.uvs.reserve(wcount * lcount * 2);
    m.colors.reserve(wcount * lcount * 4);
    m.indices.reserve(wseg * lseg * 6);

    float halfW = opts.width * 0.5f;
    float zBase = opts.stemOffset ? 0.0f : -opts.length * 0.5f;

    // Bend: the card curls forward (positive Y) along its length. We treat
    // the bend as an arc-deflection — totalAngle = bend, traced uniformly along
    // length parameter t in [0,1].
    float totalBend = opts.bend;
    bool hasBend = std::fabs(totalBend) > 1e-6f;
    float bendRadius = hasBend ? (opts.length / totalBend) : 0.0f;

    for (int j = 0; j < lcount; ++j) {
        float tj = static_cast<float>(j) / static_cast<float>(lseg);
        float zFlat = zBase + tj * opts.length;

        // Bent centerline position (in YZ plane).
        float cz, cy;
        // Tangent direction along the centerline (for normal computation).
        float tz, ty;
        if (hasBend) {
            float ang = totalBend * tj;
            cz = zBase + bendRadius * std::sin(ang);
            cy = bendRadius * (1.0f - std::cos(ang));
            tz = std::cos(ang);
            ty = std::sin(ang);
        } else {
            cz = zFlat;
            cy = 0.0f;
            tz = 1.0f;
            ty = 0.0f;
        }

        // Curl: roll around the local length axis (the tangent).
        float roll = opts.curl * tj;
        float cr = std::cos(roll);
        float sr = std::sin(roll);

        for (int i = 0; i < wcount; ++i) {
            float ti = static_cast<float>(i) / static_cast<float>(wseg);
            float xLocal = (ti - 0.5f) * opts.width;
            (void)halfW;

            // Width direction is local +X before curl. Curl rotates around the
            // tangent (cz,cy) axis. The tangent is (0, ty, tz). Local frame:
            //   tangent T = (0, ty, tz)
            //   binormal B (initially +X) = (1, 0, 0)
            //   normal N = T x B = (0*0 - tz*0, tz*1 - 0*0, 0*0 - ty*1)
            //                    = (0, tz, -ty)
            // After curl roll around T:
            //   B' = cos(roll)*B + sin(roll)*N
            //   N' = -sin(roll)*B + cos(roll)*N
            float Bx = cr * 1.0f + sr * 0.0f;
            float By = cr * 0.0f + sr * tz;
            float Bz = cr * 0.0f + sr * (-ty);

            float Nx = -sr * 1.0f + cr * 0.0f;
            float Ny = -sr * 0.0f + cr * tz;
            float Nz = -sr * 0.0f + cr * (-ty);

            float px = xLocal * Bx;
            float py = cy + xLocal * By;
            float pz = cz + xLocal * Bz;

            m.positions.push_back(px);
            m.positions.push_back(py);
            m.positions.push_back(pz);

            // Normalize normal (already unit by construction).
            float nLen = std::sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
            if (nLen > 1e-8f) { Nx /= nLen; Ny /= nLen; Nz /= nLen; }
            m.normals.push_back(Nx);
            m.normals.push_back(Ny);
            m.normals.push_back(Nz);

            // UV in atlas cell. ti -> u, tj -> v (flip v so leaf base is bottom).
            float u = u0 + ti * (u1 - u0);
            float v = v1 - tj * (v1 - v0);
            m.uvs.push_back(u);
            m.uvs.push_back(v);

            // Wind bend in R channel: 0 at base, 1 at tip. Modulate by width
            // so the card edge sways slightly more than the spine.
            float wb = tj;
            float edgeBoost = 1.0f + 0.15f * std::fabs(ti - 0.5f) * 2.0f;
            wb = std::min(1.0f, wb * edgeBoost);
            m.colors.push_back(wb);
            m.colors.push_back(0.0f);
            m.colors.push_back(0.0f);
            m.colors.push_back(1.0f);
        }
    }

    for (int j = 0; j < lseg; ++j) {
        for (int i = 0; i < wseg; ++i) {
            uint32_t a = static_cast<uint32_t>(j * wcount + i);
            uint32_t b = a + 1;
            uint32_t c = a + wcount;
            uint32_t d = c + 1;
            m.indices.push_back(a);
            m.indices.push_back(c);
            m.indices.push_back(b);
            m.indices.push_back(b);
            m.indices.push_back(c);
            m.indices.push_back(d);
        }
    }

    return m;
}

MeshData flower(const FlowerOptions& opts) {
    int petalCount = std::max(1, opts.petalCount);
    int layers = std::max(1, opts.layers);

    std::vector<MeshData> parts;
    parts.reserve(static_cast<size_t>(petalCount) * layers + 1);

    // Center dome — flattened sphere.
    {
        MeshData dome = sphere(opts.centerRadius, 12, 8);
        // Squash Y to centerHeight.
        float scaleY = (opts.centerRadius > 1e-8f)
            ? (opts.centerHeight / opts.centerRadius)
            : 1.0f;
        for (size_t i = 0; i < dome.positions.size(); i += 3) {
            dome.positions[i+1] *= scaleY;
            // Lift so base of dome sits at y=0.
            dome.positions[i+1] += opts.centerHeight * 0.5f;
        }
        // Color: paint solid centerColor into vertex colors so the merge
        // knows about the colors stream.
        dome.colors.clear();
        dome.colors.reserve(dome.vertexCount() * 4);
        for (size_t i = 0; i < dome.vertexCount(); ++i) {
            // Dome doesn't sway — windBend = 0.
            dome.colors.push_back(0.0f);
            dome.colors.push_back(0.0f);
            dome.colors.push_back(0.0f);
            dome.colors.push_back(1.0f);
        }
        // Dome has no UVs from the sphere primitive in some cases — ensure
        // it has a uv stream so merge keeps UVs across all parts.
        if (!dome.hasUVs()) {
            dome.uvs.assign(dome.vertexCount() * 2, 0.0f);
        }
        parts.push_back(std::move(dome));
    }

    for (int layer = 0; layer < layers; ++layer) {
        float layerT = (layers == 1)
            ? 0.0f
            : static_cast<float>(layer) / static_cast<float>(layers - 1);
        // Inner layers smaller, outer (layer 0) full-size.
        float scale = 1.0f - 0.4f * layerT;
        float twist = opts.layerTwist * layer;
        float yLift = opts.centerHeight * (0.4f + 0.4f * layerT);

        LeafCardOptions lo;
        lo.width = opts.petalWidth * scale;
        lo.length = opts.petalLength * scale;
        lo.bend = opts.petalBend;
        lo.curl = opts.petalCurl;
        lo.stemOffset = true;
        lo.widthSegments = 4;
        lo.lengthSegments = 8;

        for (int p = 0; p < petalCount; ++p) {
            float angle = twist + 2.0f * 3.14159265358979323846f
                * static_cast<float>(p) / static_cast<float>(petalCount);
            MeshData petal = leafCard(opts.petalShape, lo);
            // Rotate so the petal radiates outward in XZ around the Y axis,
            // pointed slightly upward via a small pre-tilt around X.
            // Pre-tilt: rotate around X so the petal lifts off the ground.
            float tilt = -0.25f - 0.15f * (1.0f - layerT);
            // Build composite rotation: first tilt around X (-pitch), then
            // yaw around Y by angle.
            // Apply analytically.
            float ct = std::cos(tilt), st = std::sin(tilt);
            for (size_t i = 0; i < petal.positions.size(); i += 3) {
                float y = petal.positions[i+1];
                float z = petal.positions[i+2];
                petal.positions[i+1] = y * ct - z * st;
                petal.positions[i+2] = y * st + z * ct;
                if (petal.hasNormals()) {
                    float ny = petal.normals[i+1];
                    float nz = petal.normals[i+2];
                    petal.normals[i+1] = ny * ct - nz * st;
                    petal.normals[i+2] = ny * st + nz * ct;
                }
            }
            float ca = std::cos(angle), sa = std::sin(angle);
            for (size_t i = 0; i < petal.positions.size(); i += 3) {
                float x = petal.positions[i];
                float z = petal.positions[i+2];
                petal.positions[i]   = x * ca + z * sa;
                petal.positions[i+2] = -x * sa + z * ca;
                petal.positions[i+1] += yLift;
                if (petal.hasNormals()) {
                    float nx = petal.normals[i];
                    float nz = petal.normals[i+2];
                    petal.normals[i]   = nx * ca + nz * sa;
                    petal.normals[i+2] = -nx * sa + nz * ca;
                }
            }
            parts.push_back(std::move(petal));
        }
    }

    return mergeMeshes(parts);
}

MeshData bladeStrip(const std::vector<Vec3>& path,
                    const BladeStripOptions& opts) {
    if (path.size() < 2) return {};
    // 4-vertex diamond in profile XY: width along X, thickness along Y. A
    // zero thickness collapses the diamond's Y vertices to the centerline,
    // which `sweep` handles fine — it just produces a degenerate ribbon
    // with two coincident edges, which the consumer can render single-sided.
    const std::vector<Vec2> profile = {
        { opts.width,        0.0f },
        { 0.0f,              opts.thickness },
        { -opts.width,       0.0f },
        { 0.0f,              -opts.thickness },
    };
    SweepOptions sopts;
    sopts.closeProfile = true;
    sopts.capStart     = opts.capStart;
    sopts.capEnd       = opts.capEnd;
    sopts.miterJoints  = opts.miterJoints;
    sopts.profileScale = opts.profileScale;
    sopts.twist        = opts.twist;
    return sweep(profile, path, sopts);
}

std::vector<Vec3> bladePath(const BladePathOptions& opts) {
    int segs = std::max(1, opts.segments);
    Vec3 dir = vnormOr(opts.tipDir, Vec3{0.0f, 1.0f, 0.0f});
    Vec3 tip = opts.base + dir * opts.length;

    // Lateral axis: project world +X perpendicular to dir; fall back to +Z
    // if dir is colinear with +X. Both candidates are then renormalized.
    Vec3 lateral;
    {
        Vec3 ax{1.0f, 0.0f, 0.0f};
        Vec3 candidate = ax - dir * vdot(ax, dir);
        if (vdot(candidate, candidate) < 1e-8f) {
            ax = Vec3{0.0f, 0.0f, 1.0f};
            candidate = ax - dir * vdot(ax, dir);
        }
        lateral = vnormOr(candidate, Vec3{1.0f, 0.0f, 0.0f});
    }

    // Quadratic Bezier control point: midpoint between base and tip, then
    // offset by `bend` on the lateral axis and `lift` along world +Y.
    Vec3 mid = opts.base + (tip - opts.base) * 0.5f;
    Vec3 ctrl = mid + lateral * opts.bend + Vec3{0.0f, opts.lift, 0.0f};

    std::vector<Vec3> out;
    out.reserve((size_t)segs + 1);
    for (int i = 0; i <= segs; ++i) {
        float t = (float)i / (float)segs;
        float u = 1.0f - t;
        // B(t) = (1-t)^2 P0 + 2(1-t)t P1 + t^2 P2
        Vec3 p = opts.base * (u * u) + ctrl * (2.0f * u * t) + tip * (t * t);
        out.push_back(p);
    }
    return out;
}

} // namespace bromesh
