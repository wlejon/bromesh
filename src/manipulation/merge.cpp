#include "bromesh/manipulation/merge.h"

namespace bromesh {

MeshData mergeMeshes(const MeshData* meshes, size_t count) {
    if (!meshes || count == 0) return {};
    if (count == 1) return meshes[0];

    // Determine which attribute streams are common to ALL inputs
    bool allHaveNormals = true;
    bool allHaveUVs = true;
    bool allHaveColors = true;
    size_t totalVerts = 0;
    size_t totalIndices = 0;

    for (size_t i = 0; i < count; ++i) {
        totalVerts += meshes[i].vertexCount();
        totalIndices += meshes[i].indices.size();
        if (!meshes[i].hasNormals()) allHaveNormals = false;
        if (!meshes[i].hasUVs()) allHaveUVs = false;
        if (!meshes[i].hasColors()) allHaveColors = false;
    }

    MeshData result;
    result.positions.reserve(totalVerts * 3);
    if (allHaveNormals) result.normals.reserve(totalVerts * 3);
    if (allHaveUVs) result.uvs.reserve(totalVerts * 2);
    if (allHaveColors) result.colors.reserve(totalVerts * 4);
    result.indices.reserve(totalIndices);

    uint32_t vertexOffset = 0;
    for (size_t i = 0; i < count; ++i) {
        const auto& m = meshes[i];

        // Append positions
        result.positions.insert(result.positions.end(),
                                m.positions.begin(), m.positions.end());

        // Append shared attribute streams
        if (allHaveNormals) {
            result.normals.insert(result.normals.end(),
                                  m.normals.begin(), m.normals.end());
        }
        if (allHaveUVs) {
            result.uvs.insert(result.uvs.end(),
                              m.uvs.begin(), m.uvs.end());
        }
        if (allHaveColors) {
            result.colors.insert(result.colors.end(),
                                 m.colors.begin(), m.colors.end());
        }

        // Append indices with offset
        for (uint32_t idx : m.indices) {
            result.indices.push_back(idx + vertexOffset);
        }

        vertexOffset += static_cast<uint32_t>(m.vertexCount());
    }

    return result;
}

MeshData mergeMeshes(const std::vector<MeshData>& meshes) {
    return mergeMeshes(meshes.data(), meshes.size());
}

} // namespace bromesh
