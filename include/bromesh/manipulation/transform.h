#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Apply a 4x4 column-major transformation matrix to positions and normals.
/// Normals are transformed by the inverse-transpose of the upper 3x3.
void transformMesh(MeshData& mesh, const float* matrix4x4);

/// Translate all vertices by (dx, dy, dz).
void translateMesh(MeshData& mesh, float dx, float dy, float dz);

/// Scale all vertices by (sx, sy, sz) about the origin.
void scaleMesh(MeshData& mesh, float sx, float sy, float sz);

/// Uniformly scale all vertices by s about the origin.
void scaleMesh(MeshData& mesh, float s);

/// Rotate all vertices and normals about an axis (ax, ay, az) by angleRadians.
/// The axis does not need to be normalized.
void rotateMesh(MeshData& mesh, float ax, float ay, float az, float angleRadians);

/// Mirror (reflect) the mesh across a plane through the origin.
/// axis: 0=YZ plane (flip X), 1=XZ plane (flip Y), 2=XY plane (flip Z).
/// Also reverses triangle winding to preserve face orientation.
void mirrorMesh(MeshData& mesh, int axis);

/// Center the mesh so its bounding-box center is at the origin.
/// Returns the translation that was applied (original center).
void centerMesh(MeshData& mesh, float* outCenter = nullptr);

} // namespace bromesh
