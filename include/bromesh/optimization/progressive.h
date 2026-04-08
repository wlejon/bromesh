#pragma once

#include "bromesh/mesh_data.h"

#include <cstdint>
#include <vector>

namespace bromesh {

/// A single edge-collapse record in a progressive mesh.
struct CollapseRecord {
    uint32_t vertexFrom;       // vertex that was removed
    uint32_t vertexTo;         // vertex it collapsed into
    float targetPosition[3];   // optimal position after collapse
    uint32_t numTrianglesRemoved; // triangles removed by this collapse (1 or 2)
};

/// Progressive mesh: stores the full-resolution mesh plus an ordered sequence
/// of edge collapses that can be replayed to reach any triangle count.
struct ProgressiveMesh {
    MeshData baseMesh;                       // full-resolution mesh (before any collapses)
    std::vector<CollapseRecord> collapses;   // ordered collapse sequence (first = coarsest removal)

    /// Number of triangles at full resolution.
    size_t maxTriangles() const { return baseMesh.triangleCount(); }

    /// Minimum triangle count reachable.
    size_t minTriangles() const;
};

/// Build a progressive mesh from the input mesh.
/// Uses quadric error metrics (QEM) to determine optimal collapse order.
/// The input mesh should have valid positions and indices.
ProgressiveMesh buildProgressiveMesh(const MeshData& mesh);

/// Extract a mesh at a specific triangle count from a progressive mesh.
/// Replays collapses from the full mesh down to the target count.
/// targetTriangles is clamped to [minTriangles, maxTriangles].
MeshData progressiveMeshAtTriangleCount(const ProgressiveMesh& pm,
                                         size_t targetTriangles);

/// Extract a mesh at a given LOD ratio (0.0 = minimum, 1.0 = full resolution).
MeshData progressiveMeshAtRatio(const ProgressiveMesh& pm, float ratio);

/// Serialize a progressive mesh to a compact binary format for streaming.
/// Format: [header][base mesh data][collapse records]
/// This enables progressive transmission: send the coarsest mesh first,
/// then stream collapse records to refine.
std::vector<uint8_t> serializeProgressiveMesh(const ProgressiveMesh& pm);

/// Deserialize a progressive mesh from binary data.
/// Returns an empty ProgressiveMesh if the data is invalid.
ProgressiveMesh deserializeProgressiveMesh(const uint8_t* data, size_t size);

} // namespace bromesh
