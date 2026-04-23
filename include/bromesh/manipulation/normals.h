#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Compute smooth (area-weighted) vertex normals from positions + indices.
/// Replaces any existing normals on the mesh.
void computeNormals(MeshData& mesh);

/// Compute flat normals (each triangle gets a uniform normal).
/// This duplicates vertices so each triangle has its own set.
MeshData computeFlatNormals(const MeshData& mesh);

/// Compute smooth normals with crease splitting. Adjacent faces whose normals
/// differ by more than angleThresholdDeg get separate vertex normals.
/// This produces smooth shading on curved surfaces and sharp edges at creases.
MeshData computeCreaseNormals(const MeshData& mesh, float angleThresholdDeg = 30.0f);

/// Compute tangent vectors (for normal mapping). Requires UVs.
/// Returns tangents as a separate array: 4 floats per vertex (xyz + handedness w).
std::vector<float> computeTangents(const MeshData& mesh);

/// Compute tangents in-place and store them on mesh.tangents (stride 4).
/// No-op if the mesh lacks UVs or normals.
void generateTangents(MeshData& mesh);

} // namespace bromesh
