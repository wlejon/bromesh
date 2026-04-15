#pragma once

#include "bromesh/mesh_data.h"
#include "bromesh/analysis/bake_texture.h"

namespace bromesh {

/// Bake a tangent-space normal map from a high-poly `reference` mesh onto
/// the UV layout of `lowPoly`.
///
/// For each texel in lowPoly's UV space, computes the world-space position
/// and tangent-space basis from the triangle that covers it, casts a ray
/// along ±normal up to `searchDistance`, samples the interpolated reference
/// normal at the hit, and packs (x,y,z) into RGB as (n*0.5+0.5). Alpha=1 on
/// covered texels, 0 on gaps.
///
/// `lowPoly` must have UVs, normals, and indices. The tangent basis is
/// derived per-triangle from UV derivatives, so no precomputed tangent
/// attribute is required.
TextureBuffer bakeNormalsFromReference(const MeshData& lowPoly,
                                        const MeshData& reference,
                                        int texWidth, int texHeight,
                                        float searchDistance = 0.0f);

/// Bake ambient occlusion against a `reference` (typically high-poly) mesh
/// onto the UV layout of `lowPoly`. Useful for transferring detail AO from
/// a dense scan onto a clean retopo cage. Returns a single-channel texture.
TextureBuffer bakeAOFromReference(const MeshData& lowPoly,
                                   const MeshData& reference,
                                   int texWidth, int texHeight,
                                   int numRays = 64,
                                   float maxDistance = 0.0f);

} // namespace bromesh
