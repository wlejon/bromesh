#include "bromesh/analysis/bake_transfer.h"
#include "bromesh/analysis/bvh.h"
#include "bromesh/analysis/raycast.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <vector>

namespace bromesh {

// ---- Barycentric coordinates in 2D (UV space) ----
static bool bary2D(float px, float py,
                    float ax, float ay, float bx, float by, float cx, float cy,
                    float& u, float& v, float& w) {
    float v0x = bx - ax, v0y = by - ay;
    float v1x = cx - ax, v1y = cy - ay;
    float v2x = px - ax, v2y = py - ay;
    float d00 = v0x*v0x + v0y*v0y;
    float d01 = v0x*v1x + v0y*v1y;
    float d11 = v1x*v1x + v1y*v1y;
    float d20 = v2x*v0x + v2y*v0y;
    float d21 = v2x*v1x + v2y*v1y;
    float denom = d00*d11 - d01*d01;
    if (std::fabs(denom) < 1e-12f) return false;
    float inv = 1.0f / denom;
    v = (d11*d20 - d01*d21) * inv;
    w = (d00*d21 - d01*d20) * inv;
    u = 1.0f - v - w;
    return (u >= -1e-4f && v >= -1e-4f && w >= -1e-4f);
}

static void normalize3(float* v) {
    float len = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if (len > 1e-8f) { v[0]/=len; v[1]/=len; v[2]/=len; }
}

static float meshBBoxDiag(const MeshData& m) {
    if (m.empty()) return 1.0f;
    float mn[3] = { m.positions[0], m.positions[1], m.positions[2] };
    float mx[3] = { mn[0], mn[1], mn[2] };
    for (size_t i = 0; i < m.positions.size(); i += 3) {
        for (int k = 0; k < 3; ++k) {
            mn[k] = std::min(mn[k], m.positions[i+k]);
            mx[k] = std::max(mx[k], m.positions[i+k]);
        }
    }
    float dx = mx[0]-mn[0], dy = mx[1]-mn[1], dz = mx[2]-mn[2];
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

// ---- Normal bake ----------------------------------------------------------

TextureBuffer bakeNormalsFromReference(const MeshData& lowPoly,
                                        const MeshData& reference,
                                        int texWidth, int texHeight,
                                        float searchDistance) {
    TextureBuffer buf;
    if (lowPoly.empty() || !lowPoly.hasUVs() || !lowPoly.hasNormals() || reference.empty())
        return buf;

    buf.width = texWidth;
    buf.height = texHeight;
    buf.channels = 4;
    buf.pixels.assign((size_t)texWidth * texHeight * 4, 0.0f);

    if (searchDistance <= 0.0f)
        searchDistance = meshBBoxDiag(lowPoly) * 0.05f;

    MeshBVH refBVH = MeshBVH::build(reference);
    const bool refHasNormals = reference.hasNormals();
    const size_t triCount = lowPoly.triangleCount();

    for (size_t t = 0; t < triCount; ++t) {
        uint32_t i0 = lowPoly.indices[t*3+0];
        uint32_t i1 = lowPoly.indices[t*3+1];
        uint32_t i2 = lowPoly.indices[t*3+2];

        const float* p0 = &lowPoly.positions[i0*3];
        const float* p1 = &lowPoly.positions[i1*3];
        const float* p2 = &lowPoly.positions[i2*3];
        const float* n0 = &lowPoly.normals[i0*3];
        const float* n1 = &lowPoly.normals[i1*3];
        const float* n2 = &lowPoly.normals[i2*3];
        float uv0u = lowPoly.uvs[i0*2+0], uv0v = lowPoly.uvs[i0*2+1];
        float uv1u = lowPoly.uvs[i1*2+0], uv1v = lowPoly.uvs[i1*2+1];
        float uv2u = lowPoly.uvs[i2*2+0], uv2v = lowPoly.uvs[i2*2+1];

        // Per-triangle tangent frame from UV derivatives.
        // (dP / dUV) * [[duv1.v, -duv1.u], [-duv2.v, duv2.u]] / det gives T/B
        float e1[3] = { p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2] };
        float e2[3] = { p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2] };
        float du1 = uv1u - uv0u, dv1 = uv1v - uv0v;
        float du2 = uv2u - uv0u, dv2 = uv2v - uv0v;
        float det = du1 * dv2 - du2 * dv1;
        float T[3], B[3];
        if (std::fabs(det) < 1e-12f) {
            // Degenerate UV — fall back to an arbitrary basis below
            T[0] = 1; T[1] = 0; T[2] = 0;
            B[0] = 0; B[1] = 1; B[2] = 0;
        } else {
            float inv = 1.0f / det;
            T[0] = (dv2 * e1[0] - dv1 * e2[0]) * inv;
            T[1] = (dv2 * e1[1] - dv1 * e2[1]) * inv;
            T[2] = (dv2 * e1[2] - dv1 * e2[2]) * inv;
            B[0] = (-du2 * e1[0] + du1 * e2[0]) * inv;
            B[1] = (-du2 * e1[1] + du1 * e2[1]) * inv;
            B[2] = (-du2 * e1[2] + du1 * e2[2]) * inv;
        }

        // Pixel-space UVs
        float x0 = uv0u * texWidth, y0 = uv0v * texHeight;
        float x1 = uv1u * texWidth, y1 = uv1v * texHeight;
        float x2 = uv2u * texWidth, y2 = uv2v * texHeight;

        int minX = std::max(0, (int)std::floor(std::min({x0, x1, x2})));
        int maxX = std::min(texWidth - 1, (int)std::ceil(std::max({x0, x1, x2})));
        int minY = std::max(0, (int)std::floor(std::min({y0, y1, y2})));
        int maxY = std::min(texHeight - 1, (int)std::ceil(std::max({y0, y1, y2})));

        for (int y = minY; y <= maxY; ++y) {
            for (int x = minX; x <= maxX; ++x) {
                float px = x + 0.5f;
                float py = y + 0.5f;
                float bu, bv, bw;
                if (!bary2D(px, py, x0, y0, x1, y1, x2, y2, bu, bv, bw)) continue;
                bu = std::max(0.0f, bu); bv = std::max(0.0f, bv); bw = std::max(0.0f, bw);
                float bsum = bu + bv + bw;
                if (bsum > 0) { bu /= bsum; bv /= bsum; bw /= bsum; }

                float pos[3] = {
                    bu*p0[0] + bv*p1[0] + bw*p2[0],
                    bu*p0[1] + bv*p1[1] + bw*p2[1],
                    bu*p0[2] + bv*p1[2] + bw*p2[2],
                };
                float N[3] = {
                    bu*n0[0] + bv*n1[0] + bw*n2[0],
                    bu*n0[1] + bv*n1[1] + bw*n2[1],
                    bu*n0[2] + bv*n1[2] + bw*n2[2],
                };
                normalize3(N);

                // Cast ray from outside-inward (offset origin along +N, cast -N)
                float origin[3] = {
                    pos[0] + N[0] * searchDistance,
                    pos[1] + N[1] * searchDistance,
                    pos[2] + N[2] * searchDistance,
                };
                float dir[3] = { -N[0], -N[1], -N[2] };
                RayHit hit = refBVH.raycast(reference, origin, dir, searchDistance * 2.0f);

                float refN[3] = { N[0], N[1], N[2] };
                if (hit.hit && refHasNormals) {
                    refN[0] = hit.normal[0];
                    refN[1] = hit.normal[1];
                    refN[2] = hit.normal[2];
                    normalize3(refN);
                }

                // Orthonormalize the T/B basis against N (Gram-Schmidt)
                float Tn[3] = { T[0], T[1], T[2] };
                float dotTN = Tn[0]*N[0] + Tn[1]*N[1] + Tn[2]*N[2];
                Tn[0] -= N[0]*dotTN; Tn[1] -= N[1]*dotTN; Tn[2] -= N[2]*dotTN;
                normalize3(Tn);
                float Bn[3] = {
                    N[1]*Tn[2] - N[2]*Tn[1],
                    N[2]*Tn[0] - N[0]*Tn[2],
                    N[0]*Tn[1] - N[1]*Tn[0],
                };
                // Flip bitangent to match UV handedness
                float handed = Bn[0]*B[0] + Bn[1]*B[1] + Bn[2]*B[2];
                if (handed < 0) { Bn[0]=-Bn[0]; Bn[1]=-Bn[1]; Bn[2]=-Bn[2]; }

                // Project refN into tangent space (TBN is column-wise [T B N])
                float tx = refN[0]*Tn[0] + refN[1]*Tn[1] + refN[2]*Tn[2];
                float ty = refN[0]*Bn[0] + refN[1]*Bn[1] + refN[2]*Bn[2];
                float tz = refN[0]*N[0]  + refN[1]*N[1]  + refN[2]*N[2];

                float* px_ = &buf.pixels[(y * texWidth + x) * 4];
                px_[0] = tx * 0.5f + 0.5f;
                px_[1] = ty * 0.5f + 0.5f;
                px_[2] = tz * 0.5f + 0.5f;
                px_[3] = 1.0f;
            }
        }
    }

    return buf;
}

// ---- AO bake against reference -------------------------------------------

static std::vector<std::array<float, 3>> hemisphereDirs(int numRays) {
    std::vector<std::array<float, 3>> out;
    out.reserve(numRays);
    float golden = 3.14159265f * (3.0f - std::sqrt(5.0f));
    for (int i = 0; i < numRays; ++i) {
        float t = (float)i / numRays;
        float cosTheta = 1.0f - t;
        float sinTheta = std::sqrt(std::max(0.0f, 1.0f - cosTheta*cosTheta));
        float phi = golden * i;
        out.push_back({ sinTheta*std::cos(phi), cosTheta, sinTheta*std::sin(phi) });
    }
    return out;
}

TextureBuffer bakeAOFromReference(const MeshData& lowPoly,
                                   const MeshData& reference,
                                   int texWidth, int texHeight,
                                   int numRays,
                                   float maxDistance) {
    TextureBuffer buf;
    if (lowPoly.empty() || !lowPoly.hasUVs() || !lowPoly.hasNormals() || reference.empty())
        return buf;

    buf.width = texWidth;
    buf.height = texHeight;
    buf.channels = 1;
    buf.pixels.assign((size_t)texWidth * texHeight, 1.0f);

    if (maxDistance <= 0.0f)
        maxDistance = meshBBoxDiag(reference) * 0.5f;

    MeshBVH bvh = MeshBVH::build(reference);
    auto samples = hemisphereDirs(numRays);
    const size_t triCount = lowPoly.triangleCount();

    for (size_t t = 0; t < triCount; ++t) {
        uint32_t i0 = lowPoly.indices[t*3+0];
        uint32_t i1 = lowPoly.indices[t*3+1];
        uint32_t i2 = lowPoly.indices[t*3+2];

        const float* p0 = &lowPoly.positions[i0*3];
        const float* p1 = &lowPoly.positions[i1*3];
        const float* p2 = &lowPoly.positions[i2*3];
        const float* n0 = &lowPoly.normals[i0*3];
        const float* n1 = &lowPoly.normals[i1*3];
        const float* n2 = &lowPoly.normals[i2*3];
        float uv0u = lowPoly.uvs[i0*2+0], uv0v = lowPoly.uvs[i0*2+1];
        float uv1u = lowPoly.uvs[i1*2+0], uv1v = lowPoly.uvs[i1*2+1];
        float uv2u = lowPoly.uvs[i2*2+0], uv2v = lowPoly.uvs[i2*2+1];

        float x0 = uv0u * texWidth, y0 = uv0v * texHeight;
        float x1 = uv1u * texWidth, y1 = uv1v * texHeight;
        float x2 = uv2u * texWidth, y2 = uv2v * texHeight;

        int minX = std::max(0, (int)std::floor(std::min({x0, x1, x2})));
        int maxX = std::min(texWidth - 1, (int)std::ceil(std::max({x0, x1, x2})));
        int minY = std::max(0, (int)std::floor(std::min({y0, y1, y2})));
        int maxY = std::min(texHeight - 1, (int)std::ceil(std::max({y0, y1, y2})));

        for (int y = minY; y <= maxY; ++y) {
            for (int x = minX; x <= maxX; ++x) {
                float px = x + 0.5f;
                float py = y + 0.5f;
                float bu, bv, bw;
                if (!bary2D(px, py, x0, y0, x1, y1, x2, y2, bu, bv, bw)) continue;
                bu = std::max(0.0f, bu); bv = std::max(0.0f, bv); bw = std::max(0.0f, bw);
                float bsum = bu + bv + bw;
                if (bsum > 0) { bu /= bsum; bv /= bsum; bw /= bsum; }

                float pos[3] = {
                    bu*p0[0] + bv*p1[0] + bw*p2[0],
                    bu*p0[1] + bv*p1[1] + bw*p2[1],
                    bu*p0[2] + bv*p1[2] + bw*p2[2],
                };
                float N[3] = {
                    bu*n0[0] + bv*n1[0] + bw*n2[0],
                    bu*n0[1] + bv*n1[1] + bw*n2[1],
                    bu*n0[2] + bv*n1[2] + bw*n2[2],
                };
                normalize3(N);

                // Tangent basis
                float Tn[3], Bn[3];
                {
                    float absX = std::fabs(N[0]), absY = std::fabs(N[1]);
                    float up[3] = { absX < absY ? 1.0f : 0.0f, absX < absY ? 0.0f : 1.0f, 0.0f };
                    Tn[0] = up[1]*N[2] - up[2]*N[1];
                    Tn[1] = up[2]*N[0] - up[0]*N[2];
                    Tn[2] = up[0]*N[1] - up[1]*N[0];
                    normalize3(Tn);
                    Bn[0] = N[1]*Tn[2] - N[2]*Tn[1];
                    Bn[1] = N[2]*Tn[0] - N[0]*Tn[2];
                    Bn[2] = N[0]*Tn[1] - N[1]*Tn[0];
                }

                int occluded = 0;
                // Offset origin slightly along normal to avoid self-hits on reference
                float ox = pos[0] + N[0] * 1e-4f;
                float oy = pos[1] + N[1] * 1e-4f;
                float oz = pos[2] + N[2] * 1e-4f;
                float origin[3] = { ox, oy, oz };
                for (const auto& s : samples) {
                    float dir[3] = {
                        s[0]*Tn[0] + s[1]*N[0] + s[2]*Bn[0],
                        s[0]*Tn[1] + s[1]*N[1] + s[2]*Bn[1],
                        s[0]*Tn[2] + s[1]*N[2] + s[2]*Bn[2],
                    };
                    if (bvh.raycastTest(reference, origin, dir, maxDistance))
                        ++occluded;
                }
                float ao = 1.0f - (float)occluded / (float)numRays;
                buf.pixels[y * texWidth + x] = ao;
            }
        }
    }

    return buf;
}

} // namespace bromesh
