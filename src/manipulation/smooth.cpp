#include "bromesh/manipulation/smooth.h"

#include <cmath>
#include <vector>

namespace bromesh {

// Build adjacency: for each vertex, list of unique neighbor vertex indices
static std::vector<std::vector<uint32_t>> buildAdjacency(const MeshData& mesh) {
    size_t vertCount = mesh.vertexCount();
    std::vector<std::vector<uint32_t>> adj(vertCount);

    size_t triCount = mesh.triangleCount();
    for (size_t t = 0; t < triCount; ++t) {
        uint32_t v[3] = {
            mesh.indices[t * 3 + 0],
            mesh.indices[t * 3 + 1],
            mesh.indices[t * 3 + 2]
        };
        for (int i = 0; i < 3; ++i) {
            uint32_t a = v[i], b = v[(i + 1) % 3];
            auto addUnique = [](std::vector<uint32_t>& list, uint32_t val) {
                for (uint32_t x : list) if (x == val) return;
                list.push_back(val);
            };
            addUnique(adj[a], b);
            addUnique(adj[b], a);
        }
    }
    return adj;
}

static void applyLaplacianStep(MeshData& mesh,
                               const std::vector<std::vector<uint32_t>>& adj,
                               float factor) {
    size_t vertCount = mesh.vertexCount();
    std::vector<float> newPos(vertCount * 3);

    for (size_t v = 0; v < vertCount; ++v) {
        const auto& neighbors = adj[v];
        if (neighbors.empty()) {
            for (int c = 0; c < 3; ++c)
                newPos[v * 3 + c] = mesh.positions[v * 3 + c];
            continue;
        }

        // Compute centroid of neighbors
        float cx = 0, cy = 0, cz = 0;
        for (uint32_t n : neighbors) {
            cx += mesh.positions[n * 3 + 0];
            cy += mesh.positions[n * 3 + 1];
            cz += mesh.positions[n * 3 + 2];
        }
        float invN = 1.0f / neighbors.size();
        cx *= invN; cy *= invN; cz *= invN;

        // Move toward centroid
        newPos[v * 3 + 0] = mesh.positions[v * 3 + 0] + factor * (cx - mesh.positions[v * 3 + 0]);
        newPos[v * 3 + 1] = mesh.positions[v * 3 + 1] + factor * (cy - mesh.positions[v * 3 + 1]);
        newPos[v * 3 + 2] = mesh.positions[v * 3 + 2] + factor * (cz - mesh.positions[v * 3 + 2]);
    }

    // Apply
    for (size_t v = 0; v < vertCount; ++v) {
        mesh.positions[v * 3 + 0] = newPos[v * 3 + 0];
        mesh.positions[v * 3 + 1] = newPos[v * 3 + 1];
        mesh.positions[v * 3 + 2] = newPos[v * 3 + 2];
    }
}

void smoothLaplacian(MeshData& mesh, float lambda, int iterations) {
    if (mesh.empty() || mesh.indices.empty() || iterations <= 0) return;

    auto adj = buildAdjacency(mesh);

    for (int i = 0; i < iterations; ++i)
        applyLaplacianStep(mesh, adj, lambda);

    // Recompute normals
    if (mesh.hasNormals()) {
        size_t vertCount = mesh.vertexCount();
        mesh.normals.assign(vertCount * 3, 0.0f);
        size_t triCount = mesh.triangleCount();
        for (size_t t = 0; t < triCount; ++t) {
            uint32_t i0 = mesh.indices[t*3+0];
            uint32_t i1 = mesh.indices[t*3+1];
            uint32_t i2 = mesh.indices[t*3+2];
            float e1x = mesh.positions[i1*3+0]-mesh.positions[i0*3+0];
            float e1y = mesh.positions[i1*3+1]-mesh.positions[i0*3+1];
            float e1z = mesh.positions[i1*3+2]-mesh.positions[i0*3+2];
            float e2x = mesh.positions[i2*3+0]-mesh.positions[i0*3+0];
            float e2y = mesh.positions[i2*3+1]-mesh.positions[i0*3+1];
            float e2z = mesh.positions[i2*3+2]-mesh.positions[i0*3+2];
            float nx = e1y*e2z - e1z*e2y;
            float ny = e1z*e2x - e1x*e2z;
            float nz = e1x*e2y - e1y*e2x;
            for (uint32_t idx : {i0, i1, i2}) {
                mesh.normals[idx*3+0] += nx;
                mesh.normals[idx*3+1] += ny;
                mesh.normals[idx*3+2] += nz;
            }
        }
        for (size_t v = 0; v < vertCount; ++v) {
            float* n = &mesh.normals[v*3];
            float len = std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
            if (len > 1e-8f) { n[0]/=len; n[1]/=len; n[2]/=len; }
        }
    }
}

void smoothTaubin(MeshData& mesh, float lambda, float mu, int iterations) {
    if (mesh.empty() || mesh.indices.empty() || iterations <= 0) return;

    auto adj = buildAdjacency(mesh);

    for (int i = 0; i < iterations; ++i) {
        applyLaplacianStep(mesh, adj, lambda); // Shrink
        applyLaplacianStep(mesh, adj, mu);     // Inflate (mu is negative)
    }

    // Recompute normals
    if (mesh.hasNormals()) {
        size_t vertCount = mesh.vertexCount();
        mesh.normals.assign(vertCount * 3, 0.0f);
        size_t triCount = mesh.triangleCount();
        for (size_t t = 0; t < triCount; ++t) {
            uint32_t i0 = mesh.indices[t*3+0];
            uint32_t i1 = mesh.indices[t*3+1];
            uint32_t i2 = mesh.indices[t*3+2];
            float e1x = mesh.positions[i1*3+0]-mesh.positions[i0*3+0];
            float e1y = mesh.positions[i1*3+1]-mesh.positions[i0*3+1];
            float e1z = mesh.positions[i1*3+2]-mesh.positions[i0*3+2];
            float e2x = mesh.positions[i2*3+0]-mesh.positions[i0*3+0];
            float e2y = mesh.positions[i2*3+1]-mesh.positions[i0*3+1];
            float e2z = mesh.positions[i2*3+2]-mesh.positions[i0*3+2];
            float nx = e1y*e2z - e1z*e2y;
            float ny = e1z*e2x - e1x*e2z;
            float nz = e1x*e2y - e1y*e2x;
            for (uint32_t idx : {i0, i1, i2}) {
                mesh.normals[idx*3+0] += nx;
                mesh.normals[idx*3+1] += ny;
                mesh.normals[idx*3+2] += nz;
            }
        }
        for (size_t v = 0; v < vertCount; ++v) {
            float* n = &mesh.normals[v*3];
            float len = std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
            if (len > 1e-8f) { n[0]/=len; n[1]/=len; n[2]/=len; }
        }
    }
}

} // namespace bromesh
