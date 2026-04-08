#include "bromesh/optimization/optimize.h"

#ifdef BROMESH_HAS_MESHOPTIMIZER
#include <meshoptimizer.h>
#include <vector>
#endif

namespace bromesh {

void optimizeVertexCache(MeshData& mesh) {
#ifdef BROMESH_HAS_MESHOPTIMIZER
    if (mesh.indices.empty()) return;

    std::vector<uint32_t> optimized(mesh.indices.size());
    meshopt_optimizeVertexCache(
        optimized.data(),
        mesh.indices.data(),
        mesh.indices.size(),
        mesh.vertexCount()
    );
    mesh.indices = std::move(optimized);
#else
    (void)mesh;
#endif
}

void optimizeVertexFetch(MeshData& mesh) {
#ifdef BROMESH_HAS_MESHOPTIMIZER
    if (mesh.indices.empty()) return;

    size_t vertexCount = mesh.vertexCount();
    size_t indexCount = mesh.indices.size();

    // Get the remap table
    std::vector<uint32_t> remap(vertexCount);
    size_t newVertCount = meshopt_optimizeVertexFetchRemap(
        remap.data(),
        mesh.indices.data(),
        indexCount,
        vertexCount
    );

    // Apply remap to each attribute stream
    auto remapStream = [&](std::vector<float>& stream, size_t stride) {
        if (stream.empty()) return;
        std::vector<float> tmp(newVertCount * stride);
        for (size_t v = 0; v < vertexCount; ++v) {
            uint32_t dst = remap[v];
            if (dst == static_cast<uint32_t>(~0u)) continue;
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
    for (size_t i = 0; i < indexCount; ++i) {
        mesh.indices[i] = remap[mesh.indices[i]];
    }
#else
    (void)mesh;
#endif
}

void optimizeOverdraw(MeshData& mesh, float threshold) {
#ifdef BROMESH_HAS_MESHOPTIMIZER
    if (mesh.indices.empty()) return;

    std::vector<uint32_t> optimized(mesh.indices.size());
    meshopt_optimizeOverdraw(
        optimized.data(),
        mesh.indices.data(),
        mesh.indices.size(),
        mesh.positions.data(),
        mesh.vertexCount(),
        sizeof(float) * 3,
        threshold
    );
    mesh.indices = std::move(optimized);
#else
    (void)mesh; (void)threshold;
#endif
}

} // namespace bromesh
