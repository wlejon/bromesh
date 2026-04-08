#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Generate a geodesic sphere by subdividing an icosahedron.
/// nsubdivisions: 0 = icosahedron (20 faces), 1 = 80 faces, 2 = 320 faces, etc.
MeshData geodesicSphere(float radius, int nsubdivisions = 2);

/// Platonic solids (unit size, centered at origin).
MeshData icosahedron();
MeshData dodecahedron();
MeshData octahedron();
MeshData tetrahedron();

/// Generate a cone along the Y axis.
/// slices: number of radial divisions. stacks: number of vertical divisions.
MeshData cone(float radius, float height, int slices = 16, int stacks = 4);

/// Generate a disc (filled circle) in the XZ plane.
MeshData disc(float radius, int slices = 16);

/// Generate a procedural rock shape.
/// seed: random seed. nsubdivisions: detail level (0-3).
MeshData rock(float radius, int seed = 42, int nsubdivisions = 2);

/// Generate a trefoil knot.
MeshData trefoilKnot(float radius, int slices = 64, int stacks = 16);

/// Generate a Klein bottle.
MeshData kleinBottle(int slices = 32, int stacks = 16);

} // namespace bromesh
