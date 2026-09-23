#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Transvoxel-style LOD chunk. This is NOT Lengyel's Transvoxel with
/// transition cells: no transition cells are generated. The chunk is
/// marching-cubed at its own LOD stride (every 2^lod samples), and then,
/// for each face whose neighbour is coarser (neighborLod > lod), every
/// vertex lying on that face has its two in-face coordinates rounded to the
/// neighbour's grid (a step of 2^neighborLod * cellSize). The seam closes
/// by moving vertices only, so snapped boundary triangles can come out
/// stretched or degenerate. A same-LOD or finer neighbour leaves the face
/// alone (the finer chunk does the snapping). Normals are recomputed after
/// the snap. Pick gridSize = k * 2^lod + 1 so the cells reach the far faces;
/// with any other size the cells stop short of them and those faces are not
/// snapped.
///
/// Sign convention: standard SDF — `field < isoLevel` is inside,
/// `field >= isoLevel` is outside. Normals point outward.
///
/// field: scalar values for this chunk's grid, row-major.
/// gridSize: cubic grid dimension (must be uniform), at least 2.
/// lod: LOD level of this chunk (0 = highest detail), in
///      [0, kTransvoxelMaxLod] with 2^lod <= gridSize - 1; anything else
///      yields an empty mesh.
/// neighborLods: LOD level of each neighbor [+X, -X, +Y, -Y, +Z, -Z].
///               -1 (or any negative) means no neighbor (chunk boundary =
///               world edge); a level above kTransvoxelMaxLod snaps as
///               kTransvoxelMaxLod.
inline constexpr int kTransvoxelMaxLod = 30;

MeshData transvoxel(const float* field, int gridSize, int lod,
                    const int neighborLods[6],
                    float isoLevel = 0.0f, float cellSize = 1.0f);

} // namespace bromesh
