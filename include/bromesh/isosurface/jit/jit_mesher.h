#pragma once

#include "bromesh/mesh_data.h"
#include "bromesh/isosurface/jit/sdf_node.h"
#include <bromath/aabb.h>

namespace bromesh {

/// Extract an isosurface triangle mesh from an SDF graph using Marching Cubes.
///
/// If computeGradients is true, calculates vertex normals analytically using
/// 1-pass forward-mode symbolic auto-diff of the compiled JIT SDF kernel.
/// If enableMacroBlockPruning is true, uses 8x8x8 interval bounds to skip empty/interior blocks.
MeshData marchingCubesFromSDF(
    const SdfGraph& graph,
    int dimX, int dimY, int dimZ,
    const bromath::AABB3& bounds,
    float isoLevel = 0.0f,
    bool closeBoundary = true,
    bool computeGradients = true,
    bool enableMacroBlockPruning = true
);

/// Extract an isosurface triangle mesh from an SDF graph using Surface Nets.
///
/// Computes smooth vertex normals analytically using 1-pass forward-mode symbolic
/// auto-diff of the compiled JIT SDF kernel.
MeshData surfaceNetsFromSDF(
    const SdfGraph& graph,
    int dimX, int dimY, int dimZ,
    const bromath::AABB3& bounds,
    float isoLevel = 0.0f,
    bool enableMacroBlockPruning = true
);

} // namespace bromesh
