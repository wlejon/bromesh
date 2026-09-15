#include "bromesh/manipulation/lathe.h"
#include "bromesh/manipulation/polygon.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace bromesh {

namespace {

constexpr float kTwoPi = 6.28318530717958647692f;

bromath::Vec3 evalPos(float r, float h, float theta, int axis) {
    float c = std::cos(theta);
    float s = std::sin(theta);
    if (axis == 0) {
        return { h, r * c, r * s };
    } else if (axis == 2) {
        return { r * c, r * s, h };
    } else { // axis == 1 (default)
        return { r * c, h, r * s };
    }
}

bromath::Vec3 evalNorm(float nr, float nh, float theta, int axis) {
    float c = std::cos(theta);
    float s = std::sin(theta);
    if (axis == 0) {
        return { nh, nr * c, nr * s };
    } else if (axis == 2) {
        return { nr * c, nr * s, nh };
    } else {
        return { nr * c, nh, nr * s };
    }
}

} // namespace

MeshData lathe(const std::vector<bromath::Vec2>& profile, const LatheOptions& opts) {
    MeshData out;
    if (profile.size() < 2) return out;

    int axis = (opts.axis >= 0 && opts.axis <= 2) ? opts.axis : 1;
    int segments = std::max(3, opts.segments);
    float totalAngle = opts.endAngle - opts.startAngle;
    if (std::fabs(totalAngle) < 1e-6f) return out;

    bool isFull = std::fabs(totalAngle) >= (kTwoPi - 1e-4f);

    // Sanitize profile: r >= 0, remove duplicate consecutive vertices
    std::vector<bromath::Vec2> prof;
    prof.reserve(profile.size());
    for (const auto& p : profile) {
        bromath::Vec2 sp{ std::fabs(p.x), p.y };
        if (!prof.empty() && bromath::vdist(sp, prof.back()) < 1e-6f) continue;
        prof.push_back(sp);
    }
    if (prof.size() < 2) return out;
    const size_t M = prof.size();

    // Cumulative arc length along profile
    std::vector<float> arc(M, 0.0f);
    for (size_t j = 1; j < M; ++j) {
        arc[j] = arc[j - 1] + bromath::vdist(prof[j], prof[j - 1]);
    }
    float totalArc = arc.back();
    if (totalArc < 1e-6f) totalArc = 1.0f;

    // 2D profile tangents and normals
    std::vector<bromath::Vec2> n2d(M);
    for (size_t j = 0; j < M; ++j) {
        bromath::Vec2 t;
        if (j == 0) {
            t = bromath::vnorm(prof[1] - prof[0]);
        } else if (j == M - 1) {
            t = bromath::vnorm(prof[M - 1] - prof[M - 2]);
        } else {
            t = bromath::vnorm(bromath::vnorm(prof[j] - prof[j - 1]) + bromath::vnorm(prof[j + 1] - prof[j]));
        }
        // Inward / outward candidate: (t.y, -t.x)
        bromath::Vec2 n{ t.y, -t.x };
        n2d[j] = bromath::vnormOr(n, { 1.0f, 0.0f });
    }

    // Orient normals so that nr > 0 at the maximum radius point
    size_t maxRIdx = 0;
    float maxR = prof[0].x;
    for (size_t j = 1; j < M; ++j) {
        if (prof[j].x > maxR) { maxR = prof[j].x; maxRIdx = j; }
    }
    if (n2d[maxRIdx].x < 0.0f) {
        for (size_t j = 0; j < M; ++j) {
            n2d[j] = -n2d[j];
        }
    }

    bool hasBottomPole = (prof[0].x < 1e-6f);
    bool hasTopPole = (prof[M - 1].x < 1e-6f);

    uint32_t bottomPoleIdx = 0;
    uint32_t topPoleIdx = 0;

    auto pushV = [&](bromath::Vec3 pos, bromath::Vec3 n, float u, float v) -> uint32_t {
        uint32_t idx = static_cast<uint32_t>(out.positions.size() / 3);
        out.positions.push_back(pos.x);
        out.positions.push_back(pos.y);
        out.positions.push_back(pos.z);
        out.normals.push_back(n.x);
        out.normals.push_back(n.y);
        out.normals.push_back(n.z);
        out.uvs.push_back(u);
        out.uvs.push_back(v);
        return idx;
    };

    auto addTri = [&](uint32_t i0, uint32_t i1, uint32_t i2) {
        out.indices.push_back(i0);
        out.indices.push_back(i1);
        out.indices.push_back(i2);
    };

    if (hasBottomPole) {
        bromath::Vec3 p = evalPos(0.0f, prof[0].y, 0.0f, axis);
        bromath::Vec3 n = evalNorm(n2d[0].x, n2d[0].y, 0.0f, axis);
        bottomPoleIdx = pushV(p, n, 0.5f, 0.0f);
    }
    if (hasTopPole) {
        bromath::Vec3 p = evalPos(0.0f, prof[M - 1].y, 0.0f, axis);
        bromath::Vec3 n = evalNorm(n2d[M - 1].x, n2d[M - 1].y, 0.0f, axis);
        topPoleIdx = pushV(p, n, 0.5f, 1.0f);
    }

    size_t startJ = hasBottomPole ? 1 : 0;
    size_t endJ = hasTopPole ? (M - 2) : (M - 1);

    // Grid of vertices: grid[s][j - startJ]
    const size_t numRows = endJ - startJ + 1;
    const int numCols = segments + 1;
    std::vector<std::vector<uint32_t>> grid(numCols, std::vector<uint32_t>(numRows));

    for (int s = 0; s < numCols; ++s) {
        float u = static_cast<float>(s) / static_cast<float>(segments);
        float theta = opts.startAngle + u * totalAngle;

        for (size_t row = 0; row < numRows; ++row) {
            size_t j = startJ + row;
            float v = arc[j] / totalArc;
            bromath::Vec3 pos;
            bromath::Vec3 norm;
            if (isFull && s == segments) {
                // Ensure bitwise identical position to s = 0 for watertight seam
                pos = evalPos(prof[j].x, prof[j].y, opts.startAngle, axis);
                norm = evalNorm(n2d[j].x, n2d[j].y, opts.startAngle, axis);
            } else {
                pos = evalPos(prof[j].x, prof[j].y, theta, axis);
                norm = evalNorm(n2d[j].x, n2d[j].y, theta, axis);
            }
            grid[s][row] = pushV(pos, norm, u, v);
        }
    }

    // Build side triangles
    for (int s = 0; s < segments; ++s) {
        int sNext = s + 1;

        // Bottom pole fan
        if (hasBottomPole) {
            uint32_t iB = grid[s][0];
            uint32_t iC = grid[sNext][0];
            bromath::Vec3 a{ out.positions[bottomPoleIdx * 3], out.positions[bottomPoleIdx * 3 + 1], out.positions[bottomPoleIdx * 3 + 2] };
            bromath::Vec3 b{ out.positions[iB * 3], out.positions[iB * 3 + 1], out.positions[iB * 3 + 2] };
            bromath::Vec3 c{ out.positions[iC * 3], out.positions[iC * 3 + 1], out.positions[iC * 3 + 2] };
            bromath::Vec3 fn = bromath::vcross(b - a, c - a);
            bromath::Vec3 nOut = evalNorm(n2d[1].x, n2d[1].y, opts.startAngle + (s + 0.5f) * (totalAngle / segments), axis);
            if (bromath::vdot(fn, nOut) >= 0.0f) {
                addTri(bottomPoleIdx, iB, iC);
            } else {
                addTri(bottomPoleIdx, iC, iB);
            }
        }

        // Middle quads
        for (size_t row = 0; row + 1 < numRows; ++row) {
            uint32_t i00 = grid[s][row];
            uint32_t i10 = grid[s][row + 1];
            uint32_t i11 = grid[sNext][row + 1];
            uint32_t i01 = grid[sNext][row];

            bromath::Vec3 a{ out.positions[i00 * 3], out.positions[i00 * 3 + 1], out.positions[i00 * 3 + 2] };
            bromath::Vec3 b{ out.positions[i10 * 3], out.positions[i10 * 3 + 1], out.positions[i10 * 3 + 2] };
            bromath::Vec3 c{ out.positions[i11 * 3], out.positions[i11 * 3 + 1], out.positions[i11 * 3 + 2] };
            bromath::Vec3 fn = bromath::vcross(b - a, c - a);
            size_t j = startJ + row;
            bromath::Vec3 nOut = evalNorm(n2d[j].x, n2d[j].y, opts.startAngle + (s + 0.5f) * (totalAngle / segments), axis);
            if (bromath::vdot(fn, nOut) >= 0.0f) {
                addTri(i00, i10, i11);
                addTri(i00, i11, i01);
            } else {
                addTri(i00, i11, i10);
                addTri(i00, i01, i11);
            }
        }

        // Top pole fan
        if (hasTopPole) {
            uint32_t iB = grid[s][numRows - 1];
            uint32_t iC = grid[sNext][numRows - 1];
            bromath::Vec3 a{ out.positions[topPoleIdx * 3], out.positions[topPoleIdx * 3 + 1], out.positions[topPoleIdx * 3 + 2] };
            bromath::Vec3 b{ out.positions[iB * 3], out.positions[iB * 3 + 1], out.positions[iB * 3 + 2] };
            bromath::Vec3 c{ out.positions[iC * 3], out.positions[iC * 3 + 1], out.positions[iC * 3 + 2] };
            bromath::Vec3 fn = bromath::vcross(b - a, c - a);
            bromath::Vec3 nOut = evalNorm(n2d[M - 2].x, n2d[M - 2].y, opts.startAngle + (s + 0.5f) * (totalAngle / segments), axis);
            if (bromath::vdot(fn, nOut) >= 0.0f) {
                addTri(topPoleIdx, iB, iC);
            } else {
                addTri(topPoleIdx, iC, iB);
            }
        }
    }

    // Partial revolution caps
    if (!isFull && (opts.capStart || opts.capEnd)) {
        std::vector<bromath::Vec2> capPoly;
        capPoly.reserve(M + 2);
        for (const auto& p : prof) capPoly.push_back(p);

        bool closedLoop = (bromath::vdist(prof.front(), prof.back()) < 1e-5f);
        if (!closedLoop) {
            if (prof.back().x > 1e-5f) capPoly.push_back({ 0.0f, prof.back().y });
            if (prof.front().x > 1e-5f) capPoly.push_back({ 0.0f, prof.front().y });
        }

        // Check 2D signed area (r = X, h = Y)
        double area2 = 0.0;
        const size_t cN = capPoly.size();
        for (size_t i = 0; i < cN; ++i) {
            size_t next = (i + 1) % cN;
            area2 += static_cast<double>(capPoly[i].x) * static_cast<double>(capPoly[next].y) -
                     static_cast<double>(capPoly[next].x) * static_cast<double>(capPoly[i].y);
        }
        if (area2 < 0.0) {
            std::reverse(capPoly.begin(), capPoly.end());
        }

        std::vector<float> outer;
        outer.reserve(capPoly.size() * 2);
        for (const auto& p : capPoly) {
            outer.push_back(p.x);
            outer.push_back(p.y);
        }

        MeshData capMesh2D = triangulatePolygon2D(outer);
        std::vector<uint32_t> capIndices;
        if (!capMesh2D.indices.empty()) {
            capIndices = capMesh2D.indices;
        } else {
            // Fallback fan from centroid
            bromath::Vec2 c2d{ 0.0f, 0.0f };
            for (const auto& p : capPoly) c2d += p;
            c2d = c2d * (1.0f / static_cast<float>(capPoly.size()));
            capPoly.push_back(c2d);
            uint32_t cIdx = static_cast<uint32_t>(capPoly.size() - 1);
            for (size_t i = 0; i < cN; ++i) {
                capIndices.push_back(cIdx);
                capIndices.push_back(static_cast<uint32_t>(i));
                capIndices.push_back(static_cast<uint32_t>((i + 1) % cN));
            }
        }

        float minR = capPoly[0].x, maxR_c = capPoly[0].x;
        float minH = capPoly[0].y, maxH_c = capPoly[0].y;
        for (const auto& p : capPoly) {
            minR = std::min(minR, p.x); maxR_c = std::max(maxR_c, p.x);
            minH = std::min(minH, p.y); maxH_c = std::max(maxH_c, p.y);
        }
        float spanR = std::max(1e-4f, maxR_c - minR);
        float spanH = std::max(1e-4f, maxH_c - minH);

        auto addCap = [&](float theta, bool isStart) {
            bromath::Vec3 tRev = evalPos(1.0f, 0.0f, theta + 0.001f, axis) - evalPos(1.0f, 0.0f, theta, axis);
            bromath::Vec3 capNormal = isStart ? -bromath::vnorm(tRev) : bromath::vnorm(tRev);

            uint32_t baseCapIdx = static_cast<uint32_t>(out.positions.size() / 3);
            for (const auto& p : capPoly) {
                bromath::Vec3 pos = evalPos(p.x, p.y, theta, axis);
                float u = (p.x - minR) / spanR;
                float v = (p.y - minH) / spanH;
                pushV(pos, capNormal, u, v);
            }

            for (size_t t = 0; t < capIndices.size(); t += 3) {
                uint32_t i0 = baseCapIdx + capIndices[t + 0];
                uint32_t i1 = baseCapIdx + capIndices[t + 1];
                uint32_t i2 = baseCapIdx + capIndices[t + 2];
                bromath::Vec3 a{ out.positions[i0 * 3], out.positions[i0 * 3 + 1], out.positions[i0 * 3 + 2] };
                bromath::Vec3 b{ out.positions[i1 * 3], out.positions[i1 * 3 + 1], out.positions[i1 * 3 + 2] };
                bromath::Vec3 c{ out.positions[i2 * 3], out.positions[i2 * 3 + 1], out.positions[i2 * 3 + 2] };
                bromath::Vec3 fn = bromath::vcross(b - a, c - a);
                if (bromath::vdot(fn, capNormal) >= 0.0f) {
                    addTri(i0, i1, i2);
                } else {
                    addTri(i0, i2, i1);
                }
            }
        };

        if (opts.capStart) addCap(opts.startAngle, true);
        if (opts.capEnd) addCap(opts.endAngle, false);
    }

    return out;
}

} // namespace bromesh
