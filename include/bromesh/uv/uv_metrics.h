#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Per-triangle UV distortion metrics.
struct UVDistortion {
    float stretch = 0.0f;       // L2 stretch (1.0 = no distortion)
    float areaDistortion = 0.0f; // ratio of UV area to 3D area (1.0 = ideal)
    float angleDistortion = 0.0f; // MIPS-like conformal error (0.0 = perfect)
};

/// Aggregate UV quality metrics over an entire mesh.
struct UVMetrics {
    float avgStretch = 0.0f;      // area-weighted average L2 stretch
    float maxStretch = 0.0f;      // worst-case triangle stretch
    float avgAreaDistortion = 0.0f;
    float maxAreaDistortion = 0.0f;
    float avgAngleDistortion = 0.0f;
    float maxAngleDistortion = 0.0f;
    float uvSpaceUsage = 0.0f;    // total UV area / bounding-box area (packing efficiency)
    size_t triangleCount = 0;
};

/// Compute per-triangle UV distortion. Returns one entry per triangle.
/// Requires mesh.uvs to be populated.
std::vector<UVDistortion> computeUVDistortion(const MeshData& mesh);

/// Compute aggregate UV quality metrics over the whole mesh.
UVMetrics measureUVQuality(const MeshData& mesh);

} // namespace bromesh
