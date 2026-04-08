#include "bromesh/manipulation/simplify.h"

#ifdef BROMESH_HAS_MESHOPTIMIZER
#include <meshoptimizer.h>
#include <algorithm>
#include <vector>
#include <cstring>
#endif

namespace bromesh {

MeshData simplify(const MeshData& mesh, float targetRatio, float targetError) {
#ifdef BROMESH_HAS_MESHOPTIMIZER
    if (mesh.empty() || mesh.indices.empty()) return {};

    size_t indexCount = mesh.indices.size();
    size_t vertexCount = mesh.vertexCount();

    // Compute target index count, rounded down to nearest multiple of 3
    size_t targetIndexCount = static_cast<size_t>(indexCount * targetRatio);
    targetIndexCount = (targetIndexCount / 3) * 3;
    if (targetIndexCount < 3) targetIndexCount = 3;

    // Run meshopt_simplify
    std::vector<uint32_t> newIndices(indexCount);
    size_t newIndexCount = meshopt_simplify(
        newIndices.data(),
        mesh.indices.data(),
        indexCount,
        mesh.positions.data(),
        vertexCount,
        sizeof(float) * 3,
        targetIndexCount,
        targetError,
        0,
        nullptr
    );
    newIndices.resize(newIndexCount);

    // Build a new mesh with only the referenced vertices
    // First, find which vertices are referenced and build a remap
    std::vector<uint32_t> remap(vertexCount, UINT32_MAX);
    uint32_t newVertCount = 0;
    for (size_t i = 0; i < newIndexCount; ++i) {
        uint32_t idx = newIndices[i];
        if (remap[idx] == UINT32_MAX) {
            remap[idx] = newVertCount++;
        }
    }

    MeshData result;
    result.positions.resize(newVertCount * 3);
    if (mesh.hasNormals()) result.normals.resize(newVertCount * 3);
    if (mesh.hasUVs()) result.uvs.resize(newVertCount * 2);
    if (mesh.hasColors()) result.colors.resize(newVertCount * 4);

    // Copy vertex attributes using remap
    for (size_t v = 0; v < vertexCount; ++v) {
        if (remap[v] == UINT32_MAX) continue;
        uint32_t dst = remap[v];
        result.positions[dst * 3 + 0] = mesh.positions[v * 3 + 0];
        result.positions[dst * 3 + 1] = mesh.positions[v * 3 + 1];
        result.positions[dst * 3 + 2] = mesh.positions[v * 3 + 2];
        if (mesh.hasNormals()) {
            result.normals[dst * 3 + 0] = mesh.normals[v * 3 + 0];
            result.normals[dst * 3 + 1] = mesh.normals[v * 3 + 1];
            result.normals[dst * 3 + 2] = mesh.normals[v * 3 + 2];
        }
        if (mesh.hasUVs()) {
            result.uvs[dst * 2 + 0] = mesh.uvs[v * 2 + 0];
            result.uvs[dst * 2 + 1] = mesh.uvs[v * 2 + 1];
        }
        if (mesh.hasColors()) {
            result.colors[dst * 4 + 0] = mesh.colors[v * 4 + 0];
            result.colors[dst * 4 + 1] = mesh.colors[v * 4 + 1];
            result.colors[dst * 4 + 2] = mesh.colors[v * 4 + 2];
            result.colors[dst * 4 + 3] = mesh.colors[v * 4 + 3];
        }
    }

    // Remap indices
    result.indices.resize(newIndexCount);
    for (size_t i = 0; i < newIndexCount; ++i) {
        result.indices[i] = remap[newIndices[i]];
    }

    return result;
#else
    (void)mesh; (void)targetRatio; (void)targetError;
    return {};
#endif
}

MeshData simplifyToTriangleCount(const MeshData& mesh, size_t targetTriangles,
                                 float targetError) {
    if (mesh.empty() || mesh.indices.empty()) return {};

    size_t currentTriangles = mesh.triangleCount();
    if (targetTriangles >= currentTriangles) return mesh;
    if (targetTriangles < 1) targetTriangles = 1;

    float ratio = static_cast<float>(targetTriangles) / static_cast<float>(currentTriangles);
    return simplify(mesh, ratio, targetError);
}

std::vector<MeshData> generateLODChain(const MeshData& mesh,
                                       const float* ratios, int count) {
#ifdef BROMESH_HAS_MESHOPTIMIZER
    std::vector<MeshData> chain;
    chain.reserve(count);
    for (int i = 0; i < count; ++i) {
        chain.push_back(simplify(mesh, ratios[i]));
    }
    return chain;
#else
    (void)mesh; (void)ratios; (void)count;
    return {};
#endif
}

} // namespace bromesh
