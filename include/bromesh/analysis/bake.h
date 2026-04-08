#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Bake ambient occlusion into vertex colors (RGBA, alpha=1).
/// Casts rays from each vertex along the hemisphere defined by its normal.
/// numRays: rays per vertex (more = better quality, slower).
/// maxDistance: maximum ray distance (0 = auto from bounding box).
/// The result is stored in mesh.colors (white=unoccluded, black=fully occluded).
void bakeAmbientOcclusion(MeshData& mesh, int numRays = 64,
                          float maxDistance = 0.0f);

/// Bake mean curvature into vertex colors.
/// Positive curvature (convex) maps toward white, negative (concave) toward black.
/// scale: multiplier on curvature before mapping to color.
void bakeCurvature(MeshData& mesh, float scale = 1.0f);

/// Bake mesh thickness into vertex colors.
/// Shoots rays inward from each vertex along the inverted normal.
/// Thicker regions are brighter. Good for subsurface scattering approximation.
void bakeThickness(MeshData& mesh, int numRays = 32,
                   float maxDistance = 0.0f);

} // namespace bromesh
