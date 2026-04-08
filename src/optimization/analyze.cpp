#include "bromesh/optimization/analyze.h"

#ifdef BROMESH_HAS_MESHOPTIMIZER
#include <meshoptimizer.h>
#endif

namespace bromesh {

VertexCacheStats analyzeVertexCache(const MeshData& mesh,
                                    unsigned int cacheSize) {
#ifdef BROMESH_HAS_MESHOPTIMIZER
    if (mesh.indices.empty()) return {};

    auto stats = meshopt_analyzeVertexCache(
        mesh.indices.data(),
        mesh.indices.size(),
        mesh.vertexCount(),
        cacheSize,
        0, 0
    );

    return {
        stats.vertices_transformed,
        stats.warps_executed,
        stats.acmr,
        stats.atvr
    };
#else
    (void)mesh; (void)cacheSize;
    return {};
#endif
}

VertexFetchStats analyzeVertexFetch(const MeshData& mesh,
                                    size_t vertexSize) {
#ifdef BROMESH_HAS_MESHOPTIMIZER
    if (mesh.indices.empty()) return {};

    auto stats = meshopt_analyzeVertexFetch(
        mesh.indices.data(),
        mesh.indices.size(),
        mesh.vertexCount(),
        vertexSize
    );

    return {
        stats.bytes_fetched,
        stats.overfetch
    };
#else
    (void)mesh; (void)vertexSize;
    return {};
#endif
}

OverdrawStats analyzeOverdraw(const MeshData& mesh) {
#ifdef BROMESH_HAS_MESHOPTIMIZER
    if (mesh.indices.empty() || mesh.positions.empty()) return {};

    auto stats = meshopt_analyzeOverdraw(
        mesh.indices.data(),
        mesh.indices.size(),
        mesh.positions.data(),
        mesh.vertexCount(),
        sizeof(float) * 3
    );

    return {
        stats.pixels_covered,
        stats.pixels_shaded,
        stats.overdraw
    };
#else
    (void)mesh;
    return {};
#endif
}

} // namespace bromesh
