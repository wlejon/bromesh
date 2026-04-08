#include "bromesh/analysis/bake_texture.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <random>
#include <vector>

namespace bromesh {

// ---- Ray-triangle intersection (Möller-Trumbore) ----

static bool rayTriangleIntersect(
    const float* orig, const float* dir,
    const float* v0, const float* v1, const float* v2,
    float& t) {

    float e1[3] = {v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]};
    float e2[3] = {v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]};

    float h[3] = {
        dir[1]*e2[2] - dir[2]*e2[1],
        dir[2]*e2[0] - dir[0]*e2[2],
        dir[0]*e2[1] - dir[1]*e2[0]
    };
    float a = e1[0]*h[0] + e1[1]*h[1] + e1[2]*h[2];
    if (std::fabs(a) < 1e-8f) return false;

    float f = 1.0f / a;
    float s[3] = {orig[0]-v0[0], orig[1]-v0[1], orig[2]-v0[2]};
    float u = f * (s[0]*h[0] + s[1]*h[1] + s[2]*h[2]);
    if (u < 0.0f || u > 1.0f) return false;

    float q[3] = {
        s[1]*e1[2] - s[2]*e1[1],
        s[2]*e1[0] - s[0]*e1[2],
        s[0]*e1[1] - s[1]*e1[0]
    };
    float v = f * (dir[0]*q[0] + dir[1]*q[1] + dir[2]*q[2]);
    if (v < 0.0f || u + v > 1.0f) return false;

    t = f * (e2[0]*q[0] + e2[1]*q[1] + e2[2]*q[2]);
    return t > 1e-5f;
}

// ---- Hemisphere sampling ----

static std::vector<std::array<float, 3>> generateHemisphereSamples(int numRays) {
    std::vector<std::array<float, 3>> dirs;
    dirs.reserve(numRays);
    float goldenAngle = 3.14159265f * (3.0f - std::sqrt(5.0f));
    for (int i = 0; i < numRays; ++i) {
        float t = static_cast<float>(i) / numRays;
        float cosTheta = 1.0f - t;
        float sinTheta = std::sqrt(1.0f - cosTheta * cosTheta);
        float phi = goldenAngle * i;
        dirs.push_back({sinTheta * std::cos(phi), cosTheta, sinTheta * std::sin(phi)});
    }
    return dirs;
}

// ---- Orthonormal basis from normal ----

static void buildBasis(const float* n, float* tangent, float* bitangent) {
    float absX = std::fabs(n[0]);
    float absY = std::fabs(n[1]);
    float up[3];
    if (absX < absY) { up[0]=1; up[1]=0; up[2]=0; }
    else             { up[0]=0; up[1]=1; up[2]=0; }

    tangent[0] = up[1]*n[2] - up[2]*n[1];
    tangent[1] = up[2]*n[0] - up[0]*n[2];
    tangent[2] = up[0]*n[1] - up[1]*n[0];
    float len = std::sqrt(tangent[0]*tangent[0]+tangent[1]*tangent[1]+tangent[2]*tangent[2]);
    if (len > 1e-8f) { tangent[0]/=len; tangent[1]/=len; tangent[2]/=len; }

    bitangent[0] = n[1]*tangent[2] - n[2]*tangent[1];
    bitangent[1] = n[2]*tangent[0] - n[0]*tangent[2];
    bitangent[2] = n[0]*tangent[1] - n[1]*tangent[0];
}

// ---- UV-space rasterization helpers ----

// Barycentric coordinates for a point (px,py) in triangle (ax,ay)-(bx,by)-(cx,cy)
// Returns true if inside, and sets u,v,w (barycentric weights for a,b,c).
static bool barycentric(float px, float py,
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

    float invDenom = 1.0f / denom;
    v = (d11*d20 - d01*d21) * invDenom;
    w = (d00*d21 - d01*d20) * invDenom;
    u = 1.0f - v - w;

    return (u >= -1e-4f && v >= -1e-4f && w >= -1e-4f);
}

// Interpolate a 3-component attribute at barycentric coords
static void interpolate3(const float* a0, const float* a1, const float* a2,
                         float u, float v, float w, float* out) {
    out[0] = u*a0[0] + v*a1[0] + w*a2[0];
    out[1] = u*a0[1] + v*a1[1] + w*a2[1];
    out[2] = u*a0[2] + v*a1[2] + w*a2[2];
}

// For each texel, find which triangle covers it and compute barycentric coords.
// Stores per-texel: triangle index (or -1), barycentric u/v/w, interpolated position, normal.
struct TexelInfo {
    int triIndex = -1;
    float baryU = 0, baryV = 0, baryW = 0;
    float pos[3] = {0, 0, 0};
    float nrm[3] = {0, 0, 0};
};

static std::vector<TexelInfo> rasterizeUVSpace(const MeshData& mesh, int w, int h) {
    std::vector<TexelInfo> texels(w * h);
    const size_t triCount = mesh.triangleCount();
    const bool hasNormals = mesh.hasNormals();

    for (size_t t = 0; t < triCount; ++t) {
        uint32_t i0 = mesh.indices[t*3+0];
        uint32_t i1 = mesh.indices[t*3+1];
        uint32_t i2 = mesh.indices[t*3+2];

        // UV coords in pixel space
        float u0 = mesh.uvs[i0*2+0] * w, v0 = mesh.uvs[i0*2+1] * h;
        float u1 = mesh.uvs[i1*2+0] * w, v1 = mesh.uvs[i1*2+1] * h;
        float u2 = mesh.uvs[i2*2+0] * w, v2 = mesh.uvs[i2*2+1] * h;

        // Bounding box in pixel coords
        int minX = std::max(0, (int)std::floor(std::min({u0, u1, u2})));
        int maxX = std::min(w - 1, (int)std::ceil(std::max({u0, u1, u2})));
        int minY = std::max(0, (int)std::floor(std::min({v0, v1, v2})));
        int maxY = std::min(h - 1, (int)std::ceil(std::max({v0, v1, v2})));

        for (int y = minY; y <= maxY; ++y) {
            for (int x = minX; x <= maxX; ++x) {
                float px = x + 0.5f;
                float py = y + 0.5f;

                float bu, bv, bw;
                if (!barycentric(px, py, u0, v0, u1, v1, u2, v2, bu, bv, bw))
                    continue;

                // Clamp for safety
                bu = std::max(0.0f, bu);
                bv = std::max(0.0f, bv);
                bw = std::max(0.0f, bw);
                float sum = bu + bv + bw;
                if (sum > 0) { bu /= sum; bv /= sum; bw /= sum; }

                auto& texel = texels[y * w + x];
                texel.triIndex = (int)t;
                texel.baryU = bu;
                texel.baryV = bv;
                texel.baryW = bw;

                interpolate3(&mesh.positions[i0*3], &mesh.positions[i1*3],
                             &mesh.positions[i2*3], bu, bv, bw, texel.pos);

                if (hasNormals) {
                    interpolate3(&mesh.normals[i0*3], &mesh.normals[i1*3],
                                 &mesh.normals[i2*3], bu, bv, bw, texel.nrm);
                    // Normalize
                    float len = std::sqrt(texel.nrm[0]*texel.nrm[0] +
                                          texel.nrm[1]*texel.nrm[1] +
                                          texel.nrm[2]*texel.nrm[2]);
                    if (len > 1e-8f) {
                        texel.nrm[0] /= len;
                        texel.nrm[1] /= len;
                        texel.nrm[2] /= len;
                    }
                }
            }
        }
    }

    return texels;
}

// ---- Auto max distance from bounding box ----

static float autoMaxDistance(const MeshData& mesh) {
    size_t vCount = mesh.vertexCount();
    float bmin[3] = {mesh.positions[0], mesh.positions[1], mesh.positions[2]};
    float bmax[3] = {mesh.positions[0], mesh.positions[1], mesh.positions[2]};
    for (size_t v = 1; v < vCount; ++v) {
        for (int c = 0; c < 3; ++c) {
            bmin[c] = std::min(bmin[c], mesh.positions[v*3+c]);
            bmax[c] = std::max(bmax[c], mesh.positions[v*3+c]);
        }
    }
    float dx = bmax[0]-bmin[0], dy = bmax[1]-bmin[1], dz = bmax[2]-bmin[2];
    return std::sqrt(dx*dx + dy*dy + dz*dz) * 0.5f;
}

// ---- Public API ----

TextureBuffer bakeAmbientOcclusionToTexture(const MeshData& mesh,
                                             int texWidth, int texHeight,
                                             int numRays, float maxDistance) {
    TextureBuffer buf;
    if (mesh.empty() || !mesh.hasUVs() || !mesh.hasNormals() || mesh.indices.empty())
        return buf;

    buf.width = texWidth;
    buf.height = texHeight;
    buf.channels = 1;
    buf.pixels.resize(texWidth * texHeight, 0.0f);

    if (maxDistance <= 0.0f) maxDistance = autoMaxDistance(mesh);

    auto texels = rasterizeUVSpace(mesh, texWidth, texHeight);
    auto samples = generateHemisphereSamples(numRays);
    const size_t triCount = mesh.triangleCount();

    for (int i = 0; i < texWidth * texHeight; ++i) {
        const auto& texel = texels[i];
        if (texel.triIndex < 0) continue;

        float tangent[3], bitangent[3];
        buildBasis(texel.nrm, tangent, bitangent);

        int hits = 0;
        for (const auto& s : samples) {
            float dir[3] = {
                s[0]*tangent[0] + s[1]*texel.nrm[0] + s[2]*bitangent[0],
                s[0]*tangent[1] + s[1]*texel.nrm[1] + s[2]*bitangent[1],
                s[0]*tangent[2] + s[1]*texel.nrm[2] + s[2]*bitangent[2]
            };
            float orig[3] = {
                texel.pos[0] + texel.nrm[0] * 1e-4f,
                texel.pos[1] + texel.nrm[1] * 1e-4f,
                texel.pos[2] + texel.nrm[2] * 1e-4f
            };

            for (size_t t = 0; t < triCount; ++t) {
                uint32_t i0 = mesh.indices[t*3+0];
                uint32_t i1 = mesh.indices[t*3+1];
                uint32_t i2 = mesh.indices[t*3+2];

                float hitT;
                if (rayTriangleIntersect(orig, dir,
                    &mesh.positions[i0*3], &mesh.positions[i1*3], &mesh.positions[i2*3],
                    hitT) && hitT < maxDistance) {
                    hits++;
                    break;
                }
            }
        }

        buf.pixels[i] = 1.0f - static_cast<float>(hits) / numRays;
    }

    return buf;
}

TextureBuffer bakeCurvatureToTexture(const MeshData& mesh,
                                      int texWidth, int texHeight,
                                      float scale) {
    TextureBuffer buf;
    if (mesh.empty() || !mesh.hasUVs() || !mesh.hasNormals() || mesh.indices.empty())
        return buf;

    // First compute per-vertex curvature
    const size_t vertCount = mesh.vertexCount();
    const size_t triCount = mesh.triangleCount();

    std::vector<std::vector<uint32_t>> adj(vertCount);
    for (size_t t = 0; t < triCount; ++t) {
        uint32_t v[3] = {mesh.indices[t*3+0], mesh.indices[t*3+1], mesh.indices[t*3+2]};
        for (int i = 0; i < 3; ++i) {
            uint32_t a = v[i], b = v[(i+1)%3];
            auto addUnique = [](std::vector<uint32_t>& list, uint32_t val) {
                for (uint32_t x : list) if (x == val) return;
                list.push_back(val);
            };
            addUnique(adj[a], b);
            addUnique(adj[b], a);
        }
    }

    std::vector<float> curvature(vertCount, 0.0f);
    float minCurv = 1e10f, maxCurv = -1e10f;

    for (size_t v = 0; v < vertCount; ++v) {
        if (adj[v].empty()) continue;
        float* pos = const_cast<float*>(&mesh.positions[v * 3]);
        float* nrm = const_cast<float*>(&mesh.normals[v * 3]);
        float meanCurv = 0.0f;
        for (uint32_t n : adj[v]) {
            float dx = mesh.positions[n*3+0] - pos[0];
            float dy = mesh.positions[n*3+1] - pos[1];
            float dz = mesh.positions[n*3+2] - pos[2];
            float dist = std::sqrt(dx*dx+dy*dy+dz*dz);
            if (dist > 1e-8f) {
                meanCurv += (dx*nrm[0] + dy*nrm[1] + dz*nrm[2]) / dist;
            }
        }
        meanCurv /= adj[v].size();
        curvature[v] = meanCurv * scale;
        minCurv = std::min(minCurv, curvature[v]);
        maxCurv = std::max(maxCurv, curvature[v]);
    }

    float range = std::max(std::fabs(minCurv), std::fabs(maxCurv));
    if (range < 1e-8f) range = 1.0f;

    // Now rasterize into texture via UV space
    buf.width = texWidth;
    buf.height = texHeight;
    buf.channels = 1;
    buf.pixels.resize(texWidth * texHeight, 0.5f); // default: flat

    auto texels = rasterizeUVSpace(mesh, texWidth, texHeight);

    for (int i = 0; i < texWidth * texHeight; ++i) {
        const auto& texel = texels[i];
        if (texel.triIndex < 0) continue;

        uint32_t i0 = mesh.indices[texel.triIndex*3+0];
        uint32_t i1 = mesh.indices[texel.triIndex*3+1];
        uint32_t i2 = mesh.indices[texel.triIndex*3+2];

        float c = texel.baryU * curvature[i0] +
                  texel.baryV * curvature[i1] +
                  texel.baryW * curvature[i2];

        float val = (c / range) * 0.5f + 0.5f;
        buf.pixels[i] = std::max(0.0f, std::min(1.0f, val));
    }

    return buf;
}

TextureBuffer bakeThicknessToTexture(const MeshData& mesh,
                                      int texWidth, int texHeight,
                                      int numRays, float maxDistance) {
    TextureBuffer buf;
    if (mesh.empty() || !mesh.hasUVs() || !mesh.hasNormals() || mesh.indices.empty())
        return buf;

    buf.width = texWidth;
    buf.height = texHeight;
    buf.channels = 1;
    buf.pixels.resize(texWidth * texHeight, 1.0f);

    if (maxDistance <= 0.0f) maxDistance = autoMaxDistance(mesh) * 2.0f;

    auto texels = rasterizeUVSpace(mesh, texWidth, texHeight);
    auto samples = generateHemisphereSamples(numRays);
    const size_t triCount = mesh.triangleCount();

    for (int i = 0; i < texWidth * texHeight; ++i) {
        const auto& texel = texels[i];
        if (texel.triIndex < 0) continue;

        float invN[3] = {-texel.nrm[0], -texel.nrm[1], -texel.nrm[2]};
        float tangent[3], bitangent[3];
        buildBasis(invN, tangent, bitangent);

        float totalThickness = 0.0f;
        int hitCount = 0;

        for (const auto& s : samples) {
            float dir[3] = {
                s[0]*tangent[0] + s[1]*invN[0] + s[2]*bitangent[0],
                s[0]*tangent[1] + s[1]*invN[1] + s[2]*bitangent[1],
                s[0]*tangent[2] + s[1]*invN[2] + s[2]*bitangent[2]
            };
            float orig[3] = {
                texel.pos[0] + invN[0] * 1e-4f,
                texel.pos[1] + invN[1] * 1e-4f,
                texel.pos[2] + invN[2] * 1e-4f
            };

            float bestT = maxDistance;
            bool hit = false;
            for (size_t t = 0; t < triCount; ++t) {
                uint32_t i0 = mesh.indices[t*3+0];
                uint32_t i1 = mesh.indices[t*3+1];
                uint32_t i2 = mesh.indices[t*3+2];

                float hitT;
                if (rayTriangleIntersect(orig, dir,
                    &mesh.positions[i0*3], &mesh.positions[i1*3], &mesh.positions[i2*3],
                    hitT) && hitT < bestT) {
                    bestT = hitT;
                    hit = true;
                }
            }
            if (hit) {
                totalThickness += bestT;
                hitCount++;
            }
        }

        buf.pixels[i] = hitCount > 0 ? totalThickness / (hitCount * maxDistance) : 1.0f;
    }

    return buf;
}

TextureBuffer bakeNormalsToTexture(const MeshData& mesh,
                                    int texWidth, int texHeight) {
    TextureBuffer buf;
    if (mesh.empty() || !mesh.hasUVs() || !mesh.hasNormals() || mesh.indices.empty())
        return buf;

    buf.width = texWidth;
    buf.height = texHeight;
    buf.channels = 4;
    buf.pixels.resize(texWidth * texHeight * 4, 0.0f);

    auto texels = rasterizeUVSpace(mesh, texWidth, texHeight);

    for (int i = 0; i < texWidth * texHeight; ++i) {
        const auto& texel = texels[i];
        if (texel.triIndex < 0) continue;

        float* px = &buf.pixels[i * 4];
        // Map normal from [-1,1] to [0,1]
        px[0] = texel.nrm[0] * 0.5f + 0.5f;
        px[1] = texel.nrm[1] * 0.5f + 0.5f;
        px[2] = texel.nrm[2] * 0.5f + 0.5f;
        px[3] = 1.0f;
    }

    return buf;
}

TextureBuffer bakePositionToTexture(const MeshData& mesh,
                                     int texWidth, int texHeight) {
    TextureBuffer buf;
    if (mesh.empty() || !mesh.hasUVs() || mesh.indices.empty())
        return buf;

    buf.width = texWidth;
    buf.height = texHeight;
    buf.channels = 4;
    buf.pixels.resize(texWidth * texHeight * 4, 0.0f);

    auto texels = rasterizeUVSpace(mesh, texWidth, texHeight);

    for (int i = 0; i < texWidth * texHeight; ++i) {
        const auto& texel = texels[i];
        if (texel.triIndex < 0) continue;

        float* px = &buf.pixels[i * 4];
        px[0] = texel.pos[0];
        px[1] = texel.pos[1];
        px[2] = texel.pos[2];
        px[3] = 1.0f;
    }

    return buf;
}

} // namespace bromesh
