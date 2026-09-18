#pragma once

#include "bromesh/mesh_data.h"
#include "bromesh/isosurface/jit/sdf_node.h"
#include <bromath/aabb.h>

namespace bromesh {

/// Extract an isosurface triangle mesh from an SDF graph using Marching Cubes.
///
/// If computeGradients is true, calculates vertex normals analytically using
/// central differences of the compiled JIT point evaluator for perfect, smooth surface normals.
MeshData marchingCubesFromSDF(
    const SdfGraph& graph,
    int dimX, int dimY, int dimZ,
    const bromath::AABB3& bounds,
    float isoLevel = 0.0f,
    bool closeBoundary = true,
    bool computeGradients = true
);

/// Extract an isosurface triangle mesh from an SDF graph using Surface Nets.
///
/// Computes smooth vertex normals analytically using central differences
/// of the compiled JIT point evaluator.
MeshData surfaceNetsFromSDF(
    const SdfGraph& graph,
    int dimX, int dimY, int dimZ,
    const bromath::AABB3& bounds,
    float isoLevel = 0.0f
);

} // namespace bromesh
