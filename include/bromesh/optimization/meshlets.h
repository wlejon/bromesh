#pragma once

#include "bromesh/mesh_data.h"
#include <vector>

namespace bromesh {

/// Culling bounds for a meshlet (bounding sphere + normal cone).
struct MeshletBounds {
    float center[3];
    float radius;
    float coneApex[3];
    float coneAxis[3];
    float coneCutoff; // cos(angle/2)
};

/// A single meshlet: a small cluster of triangles referencing a subset of vertices.
struct Meshlet {
    std::vector<uint32_t> vertices;    // indices into the original vertex buffer
    std::vector<uint8_t>  triangles;   // micro-index buffer (groups of 3, indexing into vertices[])
    MeshletBounds bounds;

    size_t triangleCount() const { return triangles.size() / 3; }
    size_t vertexCount() const { return vertices.size(); }
};

/// Parameters for meshlet generation.
struct MeshletParams {
    size_t maxVertices  = 64;   // max vertices per meshlet (<=256)
    size_t maxTriangles = 124;  // max triangles per meshlet (<=512)
    float  coneWeight   = 0.5f; // 0 = ignore cone culling, 1 = maximize cone culling
};

/// Build meshlets from a mesh. Splits geometry into GPU-friendly clusters
/// suitable for mesh shader pipelines (e.g. Nanite-style rendering).
/// Requires meshoptimizer; returns empty if unavailable.
std::vector<Meshlet> buildMeshlets(const MeshData& mesh,
                                   const MeshletParams& params = {});

} // namespace bromesh
