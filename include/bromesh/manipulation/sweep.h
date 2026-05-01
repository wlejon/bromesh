#pragma once

#include "bromesh/mesh_data.h"
#include "bromesh/procedural/vec_math.h"

#include <vector>

namespace bromesh {

struct SweepOptions {
    /// Connect last profile vertex to first to form a closed polygon.
    bool closeProfile = true;
    bool capStart = true;
    bool capEnd = true;
    /// Per-path-point scale of the profile (1.0 = original). Size 0 or 1 =
    /// constant. If size matches `path.size()`, applied per ring.
    std::vector<float> profileScale;
    /// Per-path-point twist angle in radians. Size 0 = no twist.
    std::vector<float> twist;
    /// If true, interior path vertices place the ring on the bisector plane
    /// of incoming/outgoing tangents and scale to compensate (no gaps). If
    /// false, the ring sits perpendicular to the average tangent.
    bool miterJoints = true;
};

/// Sweep a 2D profile (in the local XY plane) along a 3D `path` using
/// parallel-transport (rotation-minimizing) frames. Caps at the ends
/// triangulate the profile assuming it is convex; concave profiles will
/// produce a fan from the centroid which may overlap. Output is a
/// triangulated `MeshData` with smooth ring-averaged vertex normals.
MeshData sweep(const std::vector<Vec2>& profile,
               const std::vector<Vec3>& path,
               const SweepOptions& opts = {});

} // namespace bromesh
