#pragma once

#include "bromesh/gaussian_splat.h"
#include "bromesh/mesh_data.h"
#include <bromath/aabb.h>
#include <bromath/vec.h>
#include <vector>

namespace bromesh {

/// Apply a column-major 4x4 affine transform to a GaussianSplatCloud.
/// Transforms positions, rotates orientation quaternions, and scales linear standard deviations.
void transformSplats(GaussianSplatCloud& cloud, const float* m);

/// Translate all splats in a GaussianSplatCloud by (dx, dy, dz).
void translateSplats(GaussianSplatCloud& cloud, float dx, float dy, float dz);

/// Scale all splats in a GaussianSplatCloud by (sx, sy, sz).
void scaleSplats(GaussianSplatCloud& cloud, float sx, float sy, float sz);

struct SplatFilterOptions {
    float minOpacity = 0.005f;               // Drop splats with opacity < minOpacity
    const bromath::AABB3* cropBox = nullptr; // Optional AABB crop (drops splats outside)
    float maxScale = 0.0f;                   // If > 0, drop splats whose max linear scale > maxScale
};

/// Filter/crop splats in-place based on opacity, scale, and spatial bounding box.
void filterSplats(GaussianSplatCloud& cloud, const SplatFilterOptions& opts = {});

/// Merge multiple GaussianSplatClouds into a single unified cloud.
/// If input clouds have different SH degrees, unifies to the maximum SH degree (padding higher degrees with 0).
GaussianSplatCloud mergeSplats(const GaussianSplatCloud* clouds, size_t count);
GaussianSplatCloud mergeSplats(const std::vector<GaussianSplatCloud>& clouds);

struct MeshToSplatsOptions {
    size_t splatCount = 5000;           // Number of splat disks to sample
    float splatRadius = 0.0f;           // World radius (0 = auto-computed from surface area and splatCount)
    float opacity = 0.95f;              // Base opacity of sampled splats
    uint32_t seed = 42;                 // Deterministic random seed
};

/// Convert a MeshData with positions (and optional normals/colors) into a render-ready GaussianSplatCloud.
/// Samples points uniformly on the mesh surface, aligns splat disk orientations (rotations) with surface normals,
/// sets scales to disk radius (thin perpendicular to normal), and initializes SH DC term from vertex colors.
GaussianSplatCloud meshToSplats(const MeshData& mesh, const MeshToSplatsOptions& opts = {});

} // namespace bromesh
