#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Compute smooth (area-weighted) vertex normals from positions + indices.
/// Replaces any existing normals on the mesh.
void computeNormals(MeshData& mesh);

/// Compute flat normals (each triangle gets a uniform normal).
/// This duplicates vertices so each triangle has its own set.
MeshData computeFlatNormals(const MeshData& mesh);

/// Compute tangent vectors (for normal mapping). Requires UVs.
/// Returns tangents as a separate array: 4 floats per vertex (xyz + handedness w).
std::vector<float> computeTangents(const MeshData& mesh);

} // namespace bromesh
