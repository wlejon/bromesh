#include "bromesh/optimization/spatial.h"

#ifdef BROMESH_HAS_MESHOPTIMIZER
#include <meshoptimizer.h>
#include <vector>
#endif

namespace bromesh {

void spatialSortTriangles(MeshData& mesh) {
#ifdef BROMESH_HAS_MESHOPTIMIZER
    if (mesh.indices.empty()) return;

    std::vector<uint32_t> sorted(mesh.indices.size());
    meshopt_spatialSortTriangles(
        sorted.data(),
        mesh.indices.data(),
        mesh.indices.size(),
        mesh.positions.data(),
        mesh.vertexCount(),
        sizeof(float) * 3
    );
    mesh.indices = std::move(sorted);
#else
    (void)mesh;
#endif
}

void spatialSortVertices(MeshData& mesh) {
#ifdef BROMESH_HAS_MESHOPTIMIZER
    if (mesh.positions.empty()) return;

    size_t vertexCount = mesh.vertexCount();
    std::vector<uint32_t> remap(vertexCount);
    meshopt_spatialSortRemap(
        remap.data(),
        mesh.positions.data(),
        vertexCount,
        sizeof(float) * 3
    );

    // Apply remap to each attribute stream
    auto remapStream = [&](std::vector<float>& stream, size_t stride) {
        if (stream.empty()) return;
        std::vector<float> tmp(vertexCount * stride);
        for (size_t v = 0; v < vertexCount; ++v) {
            uint32_t dst = remap[v];
            for (size_t s = 0; s < stride; ++s) {
                tmp[dst * stride + s] = stream[v * stride + s];
            }
        }
        stream = std::move(tmp);
    };

    remapStream(mesh.positions, 3);
    remapStream(mesh.normals, 3);
    remapStream(mesh.uvs, 2);
    remapStream(mesh.colors, 4);

    // Remap indices
    for (auto& idx : mesh.indices) {
        idx = remap[idx];
    }
#else
    (void)mesh;
#endif
}

std::vector<uint32_t> generateShadowIndexBuffer(const MeshData& mesh) {
#ifdef BROMESH_HAS_MESHOPTIMIZER
    if (mesh.indices.empty()) return {};

    std::vector<uint32_t> shadow(mesh.indices.size());
    meshopt_generateShadowIndexBuffer(
        shadow.data(),
        mesh.indices.data(),
        mesh.indices.size(),
        mesh.positions.data(),
        mesh.vertexCount(),
        sizeof(float) * 3,
        sizeof(float) * 3
    );
    return shadow;
#else
    (void)mesh;
    return {};
#endif
}

} // namespace bromesh
