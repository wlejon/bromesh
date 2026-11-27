#include "bromesh/analysis/bake.h"
#include "bromesh/analysis/bvh.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <thread>
#include <vector>

namespace bromesh {

// Quasi-uniform hemisphere directions about +Y (golden-angle spiral).
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

// Orthonormal basis about n.
static void buildBasis(const float* n, float* tangent, float* bitangent) {
    float up[3];
    if (std::fabs(n[0]) < std::fabs(n[1])) { up[0] = 1; up[1] = 0; up[2] = 0; }
    else { up[0] = 0; up[1] = 1; up[2] = 0; }
    tangent[0] = up[1]*n[2] - up[2]*n[1];
    tangent[1] = up[2]*n[0] - up[0]*n[2];
    tangent[2] = up[0]*n[1] - up[1]*n[0];
    float len = std::sqrt(tangent[0]*tangent[0]+tangent[1]*tangent[1]+tangent[2]*tangent[2]);
    if (len > 1e-8f) { tangent[0]/=len; tangent[1]/=len; tangent[2]/=len; }
    bitangent[0] = n[1]*tangent[2] - n[2]*tangent[1];
    bitangent[1] = n[2]*tangent[0] - n[0]*tangent[2];
    bitangent[2] = n[0]*tangent[1] - n[1]*tangent[0];
}

static float boundsDiagonal(const MeshData& mesh) {
    size_t vertCount = mesh.vertexCount();
    float bmin[3] = {mesh.positions[0], mesh.positions[1], mesh.positions[2]};
    float bmax[3] = {bmin[0], bmin[1], bmin[2]};
    for (size_t v = 1; v < vertCount; ++v) {
        for (int c = 0; c < 3; ++c) {
            bmin[c] = std::min(bmin[c], mesh.positions[v*3+c]);
            bmax[c] = std::max(bmax[c], mesh.positions[v*3+c]);
        }
    }
    float dx = bmax[0]-bmin[0], dy = bmax[1]-bmin[1], dz = bmax[2]-bmin[2];
    return std::sqrt(dx*dx+dy*dy+dz*dz);
}

// Run fn(v) for every vertex across the hardware threads. Each call writes
// only its own vertex, so the split needs no synchronisation.
template <typename Fn>
static void forEachVertex(size_t count, Fn&& fn) {
    unsigned hw = std::max(1u, std::thread::hardware_concurrency());
    size_t workers = std::min<size_t>(hw, (count + 255) / 256);
    if (workers <= 1) { for (size_t v = 0; v < count; ++v) fn(v); return; }
    std::vector<std::thread> pool;
    pool.reserve(workers);
    for (size_t w = 0; w < workers; ++w) {
        pool.emplace_back([&, w] {
            for (size_t v = w; v < count; v += workers) fn(v);
        });
    }
    for (auto& t : pool) t.join();
}

// Hemisphere rays about `axis` (the normal, or its inverse) from each vertex,
// traced through a BVH of the mesh. The origin is lifted off the surface by
// a small fraction of the mesh's size, so a vertex never sees its own faces
// whatever the units.
void bakeAmbientOcclusion(MeshData& mesh, int numRays, float maxDistance) {
    if (mesh.empty() || mesh.indices.empty() || !mesh.hasNormals() || numRays <= 0) return;
    size_t vertCount = mesh.vertexCount();
    float diag = boundsDiagonal(mesh);
    if (maxDistance <= 0.0f) maxDistance = diag * 0.5f;
    float lift = diag * 2e-5f;
    auto samples = generateHemisphereSamples(numRays);
    MeshBVH bvh = MeshBVH::build(mesh);
    mesh.colors.resize(vertCount * 4);

    forEachVertex(vertCount, [&](size_t v) {
        const float* pos = &mesh.positions[v * 3];
        const float* nrm = &mesh.normals[v * 3];
        float tangent[3], bitangent[3];
        buildBasis(nrm, tangent, bitangent);
        float orig[3] = {pos[0] + nrm[0] * lift, pos[1] + nrm[1] * lift, pos[2] + nrm[2] * lift};
        int hits = 0;
        for (const auto& s : samples) {
            float dir[3] = {
                s[0] * tangent[0] + s[1] * nrm[0] + s[2] * bitangent[0],
                s[0] * tangent[1] + s[1] * nrm[1] + s[2] * bitangent[1],
                s[0] * tangent[2] + s[1] * nrm[2] + s[2] * bitangent[2]
            };
            if (bvh.raycastTest(mesh, orig, dir, maxDistance)) hits++;
        }
        float ao = 1.0f - static_cast<float>(hits) / numRays;
        mesh.colors[v*4+0] = ao;
        mesh.colors[v*4+1] = ao;
        mesh.colors[v*4+2] = ao;
        mesh.colors[v*4+3] = 1.0f;
    });
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
    if (mesh.empty() || mesh.indices.empty() || !mesh.hasNormals() || numRays <= 0) return;
    size_t vertCount = mesh.vertexCount();
    float diag = boundsDiagonal(mesh);
    if (maxDistance <= 0.0f) maxDistance = diag;
    float lift = diag * 2e-5f;
    auto samples = generateHemisphereSamples(numRays);
    MeshBVH bvh = MeshBVH::build(mesh);
    mesh.colors.resize(vertCount * 4);

    forEachVertex(vertCount, [&](size_t v) {
        const float* pos = &mesh.positions[v * 3];
        const float* nrm = &mesh.normals[v * 3];
        float invN[3] = {-nrm[0], -nrm[1], -nrm[2]};
        float tangent[3], bitangent[3];
        buildBasis(invN, tangent, bitangent);
        float orig[3] = {pos[0] + invN[0] * lift, pos[1] + invN[1] * lift, pos[2] + invN[2] * lift};
        float totalThickness = 0.0f;
        int hitCount = 0;
        for (const auto& s : samples) {
            float dir[3] = {
                s[0] * tangent[0] + s[1] * invN[0] + s[2] * bitangent[0],
                s[0] * tangent[1] + s[1] * invN[1] + s[2] * bitangent[1],
                s[0] * tangent[2] + s[1] * invN[2] + s[2] * bitangent[2]
            };
            RayHit hit = bvh.raycast(mesh, orig, dir, maxDistance);
            if (hit.hit) {
                totalThickness += hit.distance;
                hitCount++;
            }
        }
        float thickness = hitCount > 0 ? totalThickness / (hitCount * maxDistance) : 1.0f;
        mesh.colors[v*4+0] = thickness;
        mesh.colors[v*4+1] = thickness;
        mesh.colors[v*4+2] = thickness;
        mesh.colors[v*4+3] = 1.0f;
    });
}

} // namespace bromesh
