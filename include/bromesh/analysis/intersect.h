#pragma once

#include "bromesh/mesh_data.h"

#include <cstdint>
#include <utility>
#include <vector>

namespace bromesh {

/// A pair of triangle indices that intersect each other.
struct TrianglePair {
    uint32_t triA;
    uint32_t triB;
};

/// Find all pairs of self-intersecting triangles in a mesh.
/// Ignores triangles that share vertices (adjacent triangles).
/// Returns an empty vector if the mesh has no self-intersections.
std::vector<TrianglePair> findSelfIntersections(const MeshData& mesh);

/// Check if a mesh has any self-intersections (early-out, faster than full enumeration).
bool hasSelfIntersections(const MeshData& mesh);

/// Check if two meshes intersect each other.
/// Returns true if any triangle from mesh A intersects any triangle from mesh B.
bool meshesIntersect(const MeshData& a, const MeshData& b);

} // namespace bromesh
