#include "bromesh/manipulation/sweep.h"

#include <cmath>
#include <cstdint>

namespace bromesh {

namespace {

struct Frame {
    Vec3 t; // tangent
    Vec3 n; // normal (profile X axis in world)
    Vec3 b; // binormal (profile Y axis in world)
};

Vec3 anyPerpendicular(Vec3 t) {
    Vec3 a = (std::fabs(t.x) < 0.9f) ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
    return vnorm(vcross(t, a));
}

float scaleAt(const std::vector<float>& s, size_t i, size_t n) {
    if (s.empty()) return 1.0f;
    if (s.size() == 1) return s[0];
    if (i >= s.size()) return s.back();
    (void)n;
    return s[i];
}

float twistAt(const std::vector<float>& tw, size_t i) {
    if (tw.empty()) return 0.0f;
    if (i >= tw.size()) return tw.back();
    return tw[i];
}

} // namespace

MeshData sweep(const std::vector<Vec2>& profile,
               const std::vector<Vec3>& path,
               const SweepOptions& opts) {
    MeshData out;
    if (profile.size() < 2 || path.size() < 2) return out;

    const size_t P = profile.size();
    const size_t N = path.size();

    // Per-path-point segment tangents. tangents[i] is the direction at path[i].
    // For interior points with mitering, this is the bisector.
    std::vector<Vec3> segDir(N - 1);
    for (size_t i = 0; i + 1 < N; ++i) {
        segDir[i] = vnormOr(path[i + 1] - path[i], Vec3{0, 1, 0});
    }

    // Build frames using parallel transport.
    std::vector<Frame> frames(N);
    Vec3 t0 = segDir[0];
    Vec3 n0 = anyPerpendicular(t0);
    Vec3 b0 = vnorm(vcross(t0, n0));
    frames[0] = {t0, n0, b0};

    for (size_t i = 1; i < N; ++i) {
        Vec3 tPrev = (i == 0) ? segDir[0] : segDir[i - 1];
        Vec3 tNext = (i + 1 < N) ? segDir[i] : segDir[i - 1];
        // Rotate previous frame's basis from tPrev to tNext.
        Quat q = quatFromTo(tPrev, tNext);
        Vec3 n = vnorm(quatRotate(q, frames[i - 1].n));
        Vec3 t = (i + 1 < N) ? tNext : tPrev;
        // For miter at interior nodes, average tangents and re-orthogonalize n.
        Vec3 tFrame = t;
        if (opts.miterJoints && i + 1 < N) {
            tFrame = vnormOr(tPrev + tNext, t);
        }
        // Project n to be perpendicular to tFrame.
        n = vnormOr(n - tFrame * vdot(n, tFrame), anyPerpendicular(tFrame));
        Vec3 b = vnorm(vcross(tFrame, n));
        frames[i] = {tFrame, n, b};
    }

    // Miter scale compensation: at interior vertex with bisector, scale
    // perpendicular dimensions by 1/cos(theta/2) where theta is the angle
    // between adjacent tangents, so that contiguous segments meet without
    // gaps when the profile is extruded along the bisector plane.
    std::vector<float> miterScale(N, 1.0f);
    if (opts.miterJoints) {
        for (size_t i = 1; i + 1 < N; ++i) {
            float c = vdot(segDir[i - 1], segDir[i]);
            if (c < -0.9999f) c = -0.9999f;
            if (c > 0.9999f) c = 0.9999f;
            // cos(theta/2) = sqrt((1 + cos(theta))/2)
            float h = std::sqrt(0.5f * (1.0f + c));
            if (h > 1e-4f) miterScale[i] = 1.0f / h;
        }
    }

    // Ring vertex generation.
    out.positions.reserve(N * P * 3);
    out.normals.reserve(N * P * 3);
    out.uvs.reserve(N * P * 2);

    for (size_t i = 0; i < N; ++i) {
        const Frame& f = frames[i];
        float s = scaleAt(opts.profileScale, i, N);
        float ms = miterScale[i];
        float tw = twistAt(opts.twist, i);
        float ct = std::cos(tw);
        float st = std::sin(tw);
        for (size_t p = 0; p < P; ++p) {
            // Twist in profile XY.
            float px = profile[p].x * ct - profile[p].y * st;
            float py = profile[p].x * st + profile[p].y * ct;
            px *= s * ms;
            py *= s * ms;
            Vec3 pos = path[i] + f.n * px + f.b * py;
            out.positions.push_back(pos.x);
            out.positions.push_back(pos.y);
            out.positions.push_back(pos.z);
            // Provisional normal: outward in profile plane, will be re-summed
            // from face normals below.
            out.normals.push_back(0);
            out.normals.push_back(0);
            out.normals.push_back(0);
            float u = (P > 1) ? static_cast<float>(p) / static_cast<float>(P - 1) : 0.0f;
            float v = (N > 1) ? static_cast<float>(i) / static_cast<float>(N - 1) : 0.0f;
            out.uvs.push_back(u);
            out.uvs.push_back(v);
        }
    }

    // Triangulate side strips between consecutive rings.
    auto sideIdx = [&](size_t ring, size_t p) -> uint32_t {
        return static_cast<uint32_t>(ring * P + p);
    };
    const size_t profileEdges = opts.closeProfile ? P : (P - 1);
    out.indices.reserve((N - 1) * profileEdges * 6);
    for (size_t i = 0; i + 1 < N; ++i) {
        for (size_t p = 0; p < profileEdges; ++p) {
            size_t pNext = (p + 1) % P;
            uint32_t a = sideIdx(i, p);
            uint32_t b = sideIdx(i, pNext);
            uint32_t c = sideIdx(i + 1, p);
            uint32_t d = sideIdx(i + 1, pNext);
            out.indices.push_back(a);
            out.indices.push_back(b);
            out.indices.push_back(c);
            out.indices.push_back(b);
            out.indices.push_back(d);
            out.indices.push_back(c);
        }
    }

    auto addTri = [&](uint32_t a, uint32_t b, uint32_t c) {
        out.indices.push_back(a);
        out.indices.push_back(b);
        out.indices.push_back(c);
    };

    // Caps via centroid fans.
    if (opts.capStart || opts.capEnd) {
        // Compute profile centroid.
        Vec2 centroid{0, 0};
        for (const auto& v : profile) {
            centroid.x += v.x;
            centroid.y += v.y;
        }
        centroid.x /= static_cast<float>(P);
        centroid.y /= static_cast<float>(P);

        if (opts.capStart) {
            const Frame& f = frames.front();
            float s = scaleAt(opts.profileScale, 0, N);
            float tw = twistAt(opts.twist, 0);
            float ct = std::cos(tw), st = std::sin(tw);
            float px = (centroid.x * ct - centroid.y * st) * s;
            float py = (centroid.x * st + centroid.y * ct) * s;
            Vec3 cpos = path.front() + f.n * px + f.b * py;
            uint32_t cIdx = static_cast<uint32_t>(out.positions.size() / 3);
            out.positions.push_back(cpos.x);
            out.positions.push_back(cpos.y);
            out.positions.push_back(cpos.z);
            out.normals.push_back(-f.t.x);
            out.normals.push_back(-f.t.y);
            out.normals.push_back(-f.t.z);
            out.uvs.push_back(0.5f);
            out.uvs.push_back(0.0f);
            const size_t edges = opts.closeProfile ? P : (P - 1);
            for (size_t p = 0; p < edges; ++p) {
                size_t pNext = (p + 1) % P;
                addTri(cIdx, sideIdx(0, pNext), sideIdx(0, p));
            }
        }
        if (opts.capEnd) {
            const Frame& f = frames.back();
            float s = scaleAt(opts.profileScale, N - 1, N);
            float tw = twistAt(opts.twist, N - 1);
            float ct = std::cos(tw), st = std::sin(tw);
            float px = (centroid.x * ct - centroid.y * st) * s;
            float py = (centroid.x * st + centroid.y * ct) * s;
            Vec3 cpos = path.back() + f.n * px + f.b * py;
            uint32_t cIdx = static_cast<uint32_t>(out.positions.size() / 3);
            out.positions.push_back(cpos.x);
            out.positions.push_back(cpos.y);
            out.positions.push_back(cpos.z);
            out.normals.push_back(f.t.x);
            out.normals.push_back(f.t.y);
            out.normals.push_back(f.t.z);
            out.uvs.push_back(0.5f);
            out.uvs.push_back(1.0f);
            const size_t edges = opts.closeProfile ? P : (P - 1);
            for (size_t p = 0; p < edges; ++p) {
                size_t pNext = (p + 1) % P;
                addTri(cIdx, sideIdx(N - 1, p), sideIdx(N - 1, pNext));
            }
        }
    }

    // Compute smooth vertex normals from face normals (only ring vertices,
    // cap centroids already have an axial normal). Re-zero ring normals.
    const size_t ringVerts = N * P;
    for (size_t i = 0; i < ringVerts * 3; ++i) out.normals[i] = 0.0f;
    const size_t triCount = out.indices.size() / 3;
    for (size_t t = 0; t < triCount; ++t) {
        uint32_t i0 = out.indices[t * 3 + 0];
        uint32_t i1 = out.indices[t * 3 + 1];
        uint32_t i2 = out.indices[t * 3 + 2];
        if (i0 >= ringVerts || i1 >= ringVerts || i2 >= ringVerts) continue;
        Vec3 a{out.positions[i0*3], out.positions[i0*3+1], out.positions[i0*3+2]};
        Vec3 b{out.positions[i1*3], out.positions[i1*3+1], out.positions[i1*3+2]};
        Vec3 c{out.positions[i2*3], out.positions[i2*3+1], out.positions[i2*3+2]};
        Vec3 fn = vcross(b - a, c - a);
        for (uint32_t idx : {i0, i1, i2}) {
            out.normals[idx*3+0] += fn.x;
            out.normals[idx*3+1] += fn.y;
            out.normals[idx*3+2] += fn.z;
        }
    }
    for (size_t v = 0; v < ringVerts; ++v) {
        Vec3 n{out.normals[v*3], out.normals[v*3+1], out.normals[v*3+2]};
        Vec3 nn = vnormOr(n, {0, 1, 0});
        out.normals[v*3+0] = nn.x;
        out.normals[v*3+1] = nn.y;
        out.normals[v*3+2] = nn.z;
    }

    return out;
}

} // namespace bromesh
