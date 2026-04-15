#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

enum class ShrinkwrapMode {
    Nearest,            // move each vertex to nearest point on target
    ProjectAlongNormal, // ray-cast from vertex along its own normal (±) to target
    ProjectAlongAxis,   // ray-cast from vertex along a user axis (±) to target
};

/// Project vertices of `mesh` onto the surface of `target`.
/// `offset` pushes the projected point along the target's surface normal
/// (positive = outward; useful to keep armor slightly above the body to
/// avoid z-fighting). `maxDistance` = 0 means unlimited.
/// `axis` is only used for ProjectAlongAxis (default +Y).
/// Normals on `mesh` are NOT recomputed; call computeNormals afterward if
/// needed.
void shrinkwrap(MeshData& mesh,
                const MeshData& target,
                ShrinkwrapMode mode = ShrinkwrapMode::Nearest,
                float maxDistance = 0.0f,
                float offset = 0.0f,
                const float axis[3] = nullptr);

} // namespace bromesh
