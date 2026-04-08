#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Parameters for UV unwrapping / atlas generation.
struct UnwrapParams {
    int maxChartCount = 0;           // 0 = no limit
    float maxStretch = 0.25f;        // max allowed stretch (0-1)
    float normalDeviationWeight = 2.0f;
    float roundnessWeight = 0.01f;
    float straightnessWeight = 6.0f;
    float normalSeamWeight = 4.0f;
    float textureSeamWeight = 0.5f;
};

/// Parameters for atlas packing.
struct PackParams {
    int resolution = 0;     // 0 = auto
    int padding = 1;        // padding between charts in texels
    bool bilinear = true;   // leave space for bilinear interpolation
    bool blockAlign = false; // align charts to 4x4 blocks
    bool bruteForce = false; // brute force packing (slow but better)
};

/// Result of UV unwrapping. The mesh is modified in-place: vertices may be
/// duplicated at UV seams, and UVs are written into mesh.uvs.
/// atlasWidth/atlasHeight give the atlas dimensions if packing was performed.
struct UnwrapResult {
    int atlasWidth = 0;
    int atlasHeight = 0;
    int chartCount = 0;
    bool success = false;
};

/// Automatically unwrap a mesh and generate UV coordinates.
/// This performs chart computation, parameterization (LSCM), and atlas packing.
/// The mesh is modified: vertices at UV seams are duplicated, and mesh.uvs
/// is populated with atlas-space UV coordinates in [0,1].
/// Requires xatlas; returns unsuccessful result if unavailable.
UnwrapResult unwrapUVs(MeshData& mesh,
                       const UnwrapParams& chartParams = {},
                       const PackParams& packParams = {});

} // namespace bromesh
