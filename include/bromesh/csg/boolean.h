#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Boolean/CSG operation type.
enum class BooleanOp {
    Union,        // A + B: combine both shapes
    Difference,   // A - B: subtract B from A
    Intersection  // A ^ B: keep only overlapping volume
};

/// Perform a Boolean/CSG operation on two manifold meshes.
/// Input meshes should be closed (watertight) manifold geometry.
/// Returns the resulting mesh, or empty if the operation fails
/// (e.g. non-manifold input) or manifold library is unavailable.
MeshData booleanOp(const MeshData& a, const MeshData& b, BooleanOp op);

/// Convenience: Boolean union (A + B).
MeshData booleanUnion(const MeshData& a, const MeshData& b);

/// Convenience: Boolean difference (A - B).
MeshData booleanDifference(const MeshData& a, const MeshData& b);

/// Convenience: Boolean intersection (A ^ B).
MeshData booleanIntersection(const MeshData& a, const MeshData& b);

/// Split a mesh by a plane. Returns the two halves.
/// normal: plane normal (unit vector). offset: distance from origin along normal.
/// Both outputs may be empty if the mesh is entirely on one side.
std::pair<MeshData, MeshData> splitByPlane(const MeshData& mesh,
                                           float nx, float ny, float nz,
                                           float offset);

} // namespace bromesh
