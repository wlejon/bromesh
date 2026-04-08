#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

enum class ProjectionType {
    Box,          // Project from 6 axis directions, pick best per face
    PlanarXY,     // Project onto XY plane
    PlanarXZ,     // Project onto XZ plane
    PlanarYZ,     // Project onto YZ plane
    Cylindrical,  // Cylindrical projection around Y axis
    Spherical,    // Spherical projection from center
};

/// Generate UV coordinates by projection.
/// Replaces any existing UVs on the mesh.
/// scale: UV tiling factor (1.0 = 1 texel per world unit).
void projectUVs(MeshData& mesh, ProjectionType type, float scale = 1.0f);

} // namespace bromesh
