#include "bromesh/procedural/plants.h"

#include "bromesh/manipulation/merge.h"
#include "bromesh/manipulation/normals.h"
#include "bromesh/manipulation/sweep.h"
#include "bromesh/manipulation/transform.h"
#include "bromesh/primitives/primitives.h"
#include "bromesh/procedural/branches.h"
#include <bromath/bromath.h>

#include <algorithm>
#include <cmath>
#include <random>

namespace bromesh {

using namespace bromath;

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

// Per-shape width scale at length parameter t in [0,1]. Returns 0..1; 1 is
// the full configured width. Curves picked to match the named silhouettes:
// Oval — symmetric ellipse-like; Pointed — taper from broad base to tip;
// Petal — almond/ogive with peak ~78% length; etc.
float shapeWidthAt(LeafShape shape, float t) {
    if (t <= 0.0f) t = 0.0f;
    if (t >= 1.0f) t = 1.0f;
    switch (shape) {
    case LeafShape::Oval: {
        // sin(pi*t) — broad mid, narrow at both ends.
        return std::sin(3.14159265358979323846f * t);
    }
    case LeafShape::Pointed: {
        // Broad near base; smoothly taper to point at tip.
        // 1 at t=0.15, 0 at t=1, gentle convex curve.
        float x = std::max(0.0f, (t - 0.0f));
        return std::sqrt(1.0f - x * x);
    }
    case LeafShape::Lobed: {
        // Broad with a slight constriction near base. Approximate two
        // bumps via 0.7*sin(pi*t) + 0.3*sin(2*pi*t) clamped to [0,1].
        float a = std::sin(3.14159265358979323846f * t);
        float b = std::sin(6.28318530717958647692f * t);
        float v = 0.85f * a + 0.18f * b;
        if (v < 0.0f) v = 0.0f;
        if (v > 1.0f) v = 1.0f;
        return v;
    }
    case LeafShape::Needle: {
        // Nearly uniform thin profile, slight taper to tip.
        return 0.85f * (1.0f - 0.5f * t * t);
    }
    case LeafShape::Frond: {
        // Smooth taper from broad base to narrow tip.
        return 1.0f - 0.65f * t;
    }
    case LeafShape::Petal: {
        // Almond / ogive: narrow at base, peaks ~75% length, tapers
        // gently to a soft tip. Uses (4t(1-t))^0.6 skewed toward tip.
        float skewed = std::pow(t, 0.6f);
        float bell = 4.0f * skewed * (1.0f - skewed);
        // Lift tip a touch so it's broad-rounded, not pinched.
        float tipLift = 0.18f * std::sin(3.14159265358979323846f * t);
        float v = 0.85f * bell + tipLift;
        if (v < 0.0f) v = 0.0f;
        if (v > 1.0f) v = 1.0f;
        return v;
    }
    }
    return 1.0f;
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

        // Per-row width scale (geometric silhouette).
        float widthScale = opts.shapedSilhouette ? shapeWidthAt(shape, tj) : 1.0f;
        // Cup amount engages from base (0) to tip (full).
        float cupT = opts.cup * tj;

        for (int i = 0; i < wcount; ++i) {
            float ti = static_cast<float>(i) / static_cast<float>(wseg);
            float xUnit = (ti - 0.5f);                     // [-0.5, 0.5]
            float xLocal = xUnit * opts.width * widthScale;
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

            // Cup displacement: each width edge lifts along the local normal
            // by (xUnit)^2 * cupT * width. Symmetric across the spine so the
            // cross-section becomes a U/dish shape rather than a flat strip.
            float cupOffset = (cupT > 0.0f)
                ? cupT * xUnit * xUnit * opts.width
                : 0.0f;

            float px = xLocal * Bx + cupOffset * Nx;
            float py = cy + xLocal * By + cupOffset * Ny;
            float pz = cz + xLocal * Bz + cupOffset * Nz;

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

    // The analytical normals computed above are correct for a bent + twisted
    // strip but do not account for cup deformation or the silhouette mask.
    // Recompute when either is in play so shading matches the actual surface.
    if (opts.cup > 1e-6f || opts.shapedSilhouette) {
        computeNormals(m);
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
        float scale = 1.0f - opts.layerScaleFalloff * layerT;
        float twist = opts.layerTwist * layer;
        float yLift = opts.centerHeight *
            (opts.outerYLift + (opts.innerYLift - opts.outerYLift) * layerT);

        LeafCardOptions lo;
        lo.width = opts.petalWidth * scale;
        lo.length = opts.petalLength * scale;
        lo.bend = opts.petalBend;
        lo.curl = opts.petalCurl;
        lo.cup  = opts.petalCup;
        lo.shapedSilhouette = opts.shapedPetals;
        lo.stemOffset = true;
        lo.widthSegments = 4;
        lo.lengthSegments = 8;

        for (int p = 0; p < petalCount; ++p) {
            float angle = twist + 2.0f * 3.14159265358979323846f
                * static_cast<float>(p) / static_cast<float>(petalCount);
            MeshData petal = leafCard(opts.petalShape, lo);
            // leafCard curls the tip toward +Y for positive bend, which
            // puts its +Y normal on the INSIDE of the curl. A flower
            // viewed from outside should show the OUTSIDE (convex side)
            // of each petal, so reverse triangle winding and negate every
            // normal. The petal silhouette is unchanged; only which face
            // is "front" flips.
            for (size_t k = 0; k < petal.normals.size(); ++k) {
                petal.normals[k] = -petal.normals[k];
            }
            for (size_t t = 0; t < petal.triangleCount(); ++t) {
                std::swap(petal.indices[t * 3 + 1], petal.indices[t * 3 + 2]);
            }
            // Rotate so the petal radiates outward in XZ around the Y axis,
            // pointed slightly upward via a pre-tilt around X.
            // Lerp between outerTilt (layer 0) and innerTilt (innermost).
            float tilt = opts.outerTilt + (opts.innerTilt - opts.outerTilt) * layerT;
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

TreeResult tree(const TreeOptions& opts) {
    TreeResult out;
    if (opts.attractorCount <= 0 || opts.canopyRadius <= 0.0f) return out;

    // Uniformly sample attractors inside the canopy sphere via rejection.
    std::mt19937 rng((uint32_t)opts.seed);
    std::uniform_real_distribution<float> uni(-1.0f, 1.0f);
    std::vector<Vec3> attractors;
    attractors.reserve((size_t)opts.attractorCount);
    while ((int)attractors.size() < opts.attractorCount) {
        float x = uni(rng), y = uni(rng), z = uni(rng);
        if (x*x + y*y + z*z > 1.0f) continue;
        attractors.push_back({
            opts.canopyCenter.x + x * opts.canopyRadius,
            opts.canopyCenter.y + y * opts.canopyRadius,
            opts.canopyCenter.z + z * opts.canopyRadius,
        });
    }

    std::vector<Vec3> seeds = { opts.base };
    Vec3 initDir = vnormOr(opts.canopyCenter - opts.base, Vec3{0.0f, 1.0f, 0.0f});
    out.segments = spaceColonize(attractors, seeds, initDir, opts.colonize);
    if (out.segments.empty()) return out;

    thickenBranches(out.segments, opts.leafRadius, opts.pipeExp);
    out.branches = meshBranches(out.segments, opts.sides);
    return out;
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
