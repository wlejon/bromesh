#include "bromesh/manipulation/simplify.h"

#ifdef BROMESH_HAS_MESHOPTIMIZER
#include <meshoptimizer.h>
#include <vector>
#include <cstring>
#endif

namespace bromesh {

MeshData simplifyWithAttributes(const MeshData& mesh, float targetRatio,
                                float targetError,
                                float uvWeight, float normalWeight) {
#ifdef BROMESH_HAS_MESHOPTIMIZER
    if (mesh.empty() || mesh.indices.empty()) return {};

    size_t indexCount = mesh.indices.size();
    size_t vertexCount = mesh.vertexCount();

    size_t targetIndexCount = static_cast<size_t>(indexCount * targetRatio);
    targetIndexCount = (targetIndexCount / 3) * 3;
    if (targetIndexCount < 3) targetIndexCount = 3;

    // Build interleaved attribute array and weights
    // We include UVs (2 floats) and normals (3 floats) if present
    std::vector<float> attributes;
    std::vector<float> weights;
    size_t attrCount = 0;

    bool hasUVs = mesh.hasUVs() && uvWeight > 0.0f;
    bool hasNormals = mesh.hasNormals() && normalWeight > 0.0f;

    if (hasUVs) attrCount += 2;
    if (hasNormals) attrCount += 3;

    if (attrCount == 0) {
        // Fall back to position-only simplification
        return simplify(mesh, targetRatio, targetError);
    }

    attributes.resize(vertexCount * attrCount);
    weights.resize(attrCount);

    for (size_t v = 0; v < vertexCount; ++v) {
        size_t ai = v * attrCount;
        size_t off = 0;
        if (hasUVs) {
            attributes[ai + off + 0] = mesh.uvs[v * 2 + 0];
            attributes[ai + off + 1] = mesh.uvs[v * 2 + 1];
            off += 2;
        }
        if (hasNormals) {
            attributes[ai + off + 0] = mesh.normals[v * 3 + 0];
            attributes[ai + off + 1] = mesh.normals[v * 3 + 1];
            attributes[ai + off + 2] = mesh.normals[v * 3 + 2];
            off += 3;
        }
    }

    // Fill weights (only once, outside the vertex loop)
    {
        size_t off = 0;
        if (hasUVs) {
            weights[off + 0] = uvWeight;
            weights[off + 1] = uvWeight;
            off += 2;
        }
        if (hasNormals) {
            weights[off + 0] = normalWeight;
            weights[off + 1] = normalWeight;
            weights[off + 2] = normalWeight;
        }
    }

    std::vector<uint32_t> newIndices(indexCount);
    size_t newIndexCount = meshopt_simplifyWithAttributes(
        newIndices.data(),
        mesh.indices.data(),
        indexCount,
        mesh.positions.data(),
        vertexCount,
        sizeof(float) * 3,
        attributes.data(),
        sizeof(float) * attrCount,
        weights.data(),
        attrCount,
        nullptr, // no vertex lock
        targetIndexCount,
        targetError,
        0,
        nullptr
    );
    newIndices.resize(newIndexCount);

    // Build compacted output mesh (same pattern as simplify)
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

    result.indices.resize(newIndexCount);
    for (size_t i = 0; i < newIndexCount; ++i) {
        result.indices[i] = remap[newIndices[i]];
    }

    return result;
#else
    (void)mesh; (void)targetRatio; (void)targetError;
    (void)uvWeight; (void)normalWeight;
    return {};
#endif
}

} // namespace bromesh
