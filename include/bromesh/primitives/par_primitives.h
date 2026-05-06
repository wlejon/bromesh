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

/// Generate a cone along the Y axis. The base disc sits at Y=0 with radius
/// `radius`; the apex sits at Y=`height`. With `capBase=true` a fan cap is
/// added at Y=0 with normals pointing -Y; otherwise the bottom is open.
/// slices: number of radial divisions. stacks: number of vertical divisions.
MeshData cone(float radius, float height, int slices = 16, int stacks = 4,
              bool capBase = false);

/// Generate a disc (filled circle) in the XZ plane.
MeshData disc(float radius, int slices = 16);

/// Generate a procedural rock shape.
/// seed: random seed. nsubdivisions: detail level (0-3).
MeshData rock(float radius, int seed = 42, int nsubdivisions = 2);

/// Noise-displaced sphere with non-uniform scale and an optional translation
/// baked into the geometry — saves a `.scale().translate()` round-trip when
/// stamping many "canopy blob"-style shapes. Equivalent to:
///   rock(radius, seed, nsubdivisions)
///     .scale(scaleX, scaleY, scaleZ)
///     .translate(centerX, centerY, centerZ)
/// but normals are recomputed once at the end and the intermediate copies
/// are skipped.
MeshData blob(float radius,
              int   seed           = 42,
              int   nsubdivisions  = 2,
              float scaleX         = 1.0f,
              float scaleY         = 1.0f,
              float scaleZ         = 1.0f,
              float centerX        = 0.0f,
              float centerY        = 0.0f,
              float centerZ        = 0.0f);

/// Generate a trefoil knot.
MeshData trefoilKnot(float radius, int slices = 64, int stacks = 16);

/// Generate a Klein bottle.
MeshData kleinBottle(int slices = 32, int stacks = 16);

} // namespace bromesh
