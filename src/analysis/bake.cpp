#include "bromesh/analysis/bake.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <random>
#include <vector>

namespace bromesh {

// Simple ray-triangle intersection (Möller-Trumbore)
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

// Generate hemisphere sample directions using cosine-weighted distribution
static std::vector<std::array<float, 3>> generateHemisphereSamples(int numRays) {
    std::vector<std::array<float, 3>> dirs;
    dirs.reserve(numRays);

    // Use golden ratio spiral for quasi-uniform hemisphere sampling
    float goldenAngle = 3.14159265f * (3.0f - std::sqrt(5.0f));

    for (int i = 0; i < numRays; ++i) {
        float t = static_cast<float>(i) / numRays;
        float cosTheta = 1.0f - t; // Hemisphere: 0 to 1
        float sinTheta = std::sqrt(1.0f - cosTheta * cosTheta);
        float phi = goldenAngle * i;

        dirs.push_back({sinTheta * std::cos(phi),
                        cosTheta,
                        sinTheta * std::sin(phi)});
    }
    return dirs;
}

// Build orthonormal basis from normal
static void buildBasis(const float* n, float* tangent, float* bitangent) {
    float absX = std::fabs(n[0]);
    float absY = std::fabs(n[1]);
    float up[3];
    if (absX < absY) {
        up[0] = 1; up[1] = 0; up[2] = 0;
    } else {
        up[0] = 0; up[1] = 1; up[2] = 0;
    }

    // tangent = normalize(cross(up, n))
    tangent[0] = up[1]*n[2] - up[2]*n[1];
    tangent[1] = up[2]*n[0] - up[0]*n[2];
    tangent[2] = up[0]*n[1] - up[1]*n[0];
    float len = std::sqrt(tangent[0]*tangent[0]+tangent[1]*tangent[1]+tangent[2]*tangent[2]);
    if (len > 1e-8f) { tangent[0]/=len; tangent[1]/=len; tangent[2]/=len; }

    // bitangent = cross(n, tangent)
    bitangent[0] = n[1]*tangent[2] - n[2]*tangent[1];
    bitangent[1] = n[2]*tangent[0] - n[0]*tangent[2];
    bitangent[2] = n[0]*tangent[1] - n[1]*tangent[0];
}

void bakeAmbientOcclusion(MeshData& mesh, int numRays, float maxDistance) {
    if (mesh.empty() || mesh.indices.empty() || !mesh.hasNormals()) return;

    size_t vertCount = mesh.vertexCount();
    size_t triCount = mesh.triangleCount();

    // Auto max distance from bounding box
    if (maxDistance <= 0.0f) {
        float bmin[3] = {mesh.positions[0], mesh.positions[1], mesh.positions[2]};
        float bmax[3] = {mesh.positions[0], mesh.positions[1], mesh.positions[2]};
        for (size_t v = 1; v < vertCount; ++v) {
            for (int c = 0; c < 3; ++c) {
                bmin[c] = std::min(bmin[c], mesh.positions[v*3+c]);
                bmax[c] = std::max(bmax[c], mesh.positions[v*3+c]);
            }
        }
        float dx = bmax[0]-bmin[0], dy = bmax[1]-bmin[1], dz = bmax[2]-bmin[2];
        maxDistance = std::sqrt(dx*dx+dy*dy+dz*dz) * 0.5f;
    }

    auto samples = generateHemisphereSamples(numRays);

    mesh.colors.resize(vertCount * 4);

    for (size_t v = 0; v < vertCount; ++v) {
        float* pos = &mesh.positions[v * 3];
        float* nrm = &mesh.normals[v * 3];

        float tangent[3], bitangent[3];
        buildBasis(nrm, tangent, bitangent);

        int hits = 0;
        for (const auto& s : samples) {
            // Transform sample to world space
            float dir[3] = {
                s[0] * tangent[0] + s[1] * nrm[0] + s[2] * bitangent[0],
                s[0] * tangent[1] + s[1] * nrm[1] + s[2] * bitangent[1],
                s[0] * tangent[2] + s[1] * nrm[2] + s[2] * bitangent[2]
            };

            // Offset origin slightly along normal to avoid self-intersection
            float orig[3] = {
                pos[0] + nrm[0] * 1e-4f,
                pos[1] + nrm[1] * 1e-4f,
                pos[2] + nrm[2] * 1e-4f
            };

            // Test against all triangles
            for (size_t t = 0; t < triCount; ++t) {
                uint32_t i0 = mesh.indices[t*3+0];
                uint32_t i1 = mesh.indices[t*3+1];
                uint32_t i2 = mesh.indices[t*3+2];

                // Skip triangles containing this vertex
                if (i0 == v || i1 == v || i2 == v) continue;

                float hitT;
                if (rayTriangleIntersect(orig, dir,
                    &mesh.positions[i0*3], &mesh.positions[i1*3], &mesh.positions[i2*3],
                    hitT) && hitT < maxDistance) {
                    hits++;
                    break; // Only need first hit
                }
            }
        }

        float ao = 1.0f - static_cast<float>(hits) / numRays;
        mesh.colors[v*4+0] = ao;
        mesh.colors[v*4+1] = ao;
        mesh.colors[v*4+2] = ao;
        mesh.colors[v*4+3] = 1.0f;
    }
}

void bakeCurvature(MeshData& mesh, float scale) {
    if (mesh.empty() || mesh.indices.empty() || !mesh.hasNormals()) return;

    size_t vertCount = mesh.vertexCount();

    // Build adjacency
    std::vector<std::vector<uint32_t>> adj(vertCount);
    size_t triCount = mesh.triangleCount();
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

    // Compute mean curvature per vertex:
    // H ≈ (1/2n) * Σ dot(p_neighbor - p, normal)
    // Normalized to [-1, 1] range via scale
    std::vector<float> curvature(vertCount, 0.0f);
    float minCurv = 1e10f, maxCurv = -1e10f;

    for (size_t v = 0; v < vertCount; ++v) {
        if (adj[v].empty()) continue;

        float* pos = &mesh.positions[v * 3];
        float* nrm = &mesh.normals[v * 3];
        float meanCurv = 0.0f;

        for (uint32_t n : adj[v]) {
            float dx = mesh.positions[n*3+0] - pos[0];
            float dy = mesh.positions[n*3+1] - pos[1];
            float dz = mesh.positions[n*3+2] - pos[2];
            float dist = std::sqrt(dx*dx+dy*dy+dz*dz);
            if (dist > 1e-8f) {
                // Signed curvature: projection of edge onto normal
                float proj = (dx*nrm[0] + dy*nrm[1] + dz*nrm[2]) / dist;
                meanCurv += proj;
            }
        }
        meanCurv /= adj[v].size();
        curvature[v] = meanCurv * scale;
        minCurv = std::min(minCurv, curvature[v]);
        maxCurv = std::max(maxCurv, curvature[v]);
    }

    // Map to color: negative curvature -> black, positive -> white
    mesh.colors.resize(vertCount * 4);
    float range = std::max(std::fabs(minCurv), std::fabs(maxCurv));
    if (range < 1e-8f) range = 1.0f;

    for (size_t v = 0; v < vertCount; ++v) {
        float val = (curvature[v] / range) * 0.5f + 0.5f;
        val = std::max(0.0f, std::min(1.0f, val));
        mesh.colors[v*4+0] = val;
        mesh.colors[v*4+1] = val;
        mesh.colors[v*4+2] = val;
        mesh.colors[v*4+3] = 1.0f;
    }
}

void bakeThickness(MeshData& mesh, int numRays, float maxDistance) {
    if (mesh.empty() || mesh.indices.empty() || !mesh.hasNormals()) return;

    size_t vertCount = mesh.vertexCount();
    size_t triCount = mesh.triangleCount();

    if (maxDistance <= 0.0f) {
        float bmin[3] = {mesh.positions[0], mesh.positions[1], mesh.positions[2]};
        float bmax[3] = {mesh.positions[0], mesh.positions[1], mesh.positions[2]};
        for (size_t v = 1; v < vertCount; ++v) {
            for (int c = 0; c < 3; ++c) {
                bmin[c] = std::min(bmin[c], mesh.positions[v*3+c]);
                bmax[c] = std::max(bmax[c], mesh.positions[v*3+c]);
            }
        }
        float dx = bmax[0]-bmin[0], dy = bmax[1]-bmin[1], dz = bmax[2]-bmin[2];
        maxDistance = std::sqrt(dx*dx+dy*dy+dz*dz);
    }

    // Hemisphere samples in the inverted normal direction
    auto samples = generateHemisphereSamples(numRays);

    mesh.colors.resize(vertCount * 4);

    for (size_t v = 0; v < vertCount; ++v) {
        float* pos = &mesh.positions[v * 3];
        float* nrm = &mesh.normals[v * 3];

        // Invert normal for thickness rays
        float invN[3] = {-nrm[0], -nrm[1], -nrm[2]};
        float tangent[3], bitangent[3];
        buildBasis(invN, tangent, bitangent);

        float totalThickness = 0.0f;
        int hitCount = 0;

        for (const auto& s : samples) {
            float dir[3] = {
                s[0] * tangent[0] + s[1] * invN[0] + s[2] * bitangent[0],
                s[0] * tangent[1] + s[1] * invN[1] + s[2] * bitangent[1],
                s[0] * tangent[2] + s[1] * invN[2] + s[2] * bitangent[2]
            };

            float orig[3] = {
                pos[0] + invN[0] * 1e-4f,
                pos[1] + invN[1] * 1e-4f,
                pos[2] + invN[2] * 1e-4f
            };

            float bestT = maxDistance;
            bool hit = false;
            for (size_t t = 0; t < triCount; ++t) {
                uint32_t i0 = mesh.indices[t*3+0];
                uint32_t i1 = mesh.indices[t*3+1];
                uint32_t i2 = mesh.indices[t*3+2];
                if (i0 == v || i1 == v || i2 == v) continue;

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

        float thickness = hitCount > 0 ? totalThickness / (hitCount * maxDistance) : 1.0f;
        mesh.colors[v*4+0] = thickness;
        mesh.colors[v*4+1] = thickness;
        mesh.colors[v*4+2] = thickness;
        mesh.colors[v*4+3] = 1.0f;
    }
}

} // namespace bromesh
