#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Vertex cache statistics.
struct VertexCacheStats {
    unsigned int verticesTransformed;
    unsigned int warpsExecuted;
    float acmr;  // transformed vertices / triangle count (best ~0.5)
    float atvr;  // transformed vertices / vertex count (best 1.0)
};

/// Vertex fetch statistics.
struct VertexFetchStats {
    unsigned int bytesFetched;
    float overfetch;  // fetched bytes / vertex buffer size (best 1.0)
};

/// Overdraw statistics.
struct OverdrawStats {
    unsigned int pixelsCovered;
    unsigned int pixelsShaded;
    float overdraw;  // shaded / covered (best 1.0)
};

/// Analyze vertex cache efficiency of a mesh.
/// cache_size: simulated FIFO cache size (default 16).
VertexCacheStats analyzeVertexCache(const MeshData& mesh,
                                    unsigned int cacheSize = 16);

/// Analyze vertex fetch efficiency of a mesh.
/// vertexSize: size of each vertex in bytes for simulation.
VertexFetchStats analyzeVertexFetch(const MeshData& mesh,
                                    size_t vertexSize = 32);

/// Analyze overdraw of a mesh using software rasterization.
OverdrawStats analyzeOverdraw(const MeshData& mesh);

} // namespace bromesh
