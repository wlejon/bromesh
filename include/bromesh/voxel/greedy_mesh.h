#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Greedy meshing for discrete voxel volumes.
/// Merges coplanar adjacent faces into larger quads to minimize triangle count.
///
/// voxels: 3D grid of voxel IDs (0 = empty, nonzero = solid material).
/// gridX/Y/Z: grid dimensions.
/// cellSize: world-space size of each voxel (default 1.0).
///
/// Output includes vertex colors derived from the palette if provided.
/// palette: RGBA float array, 4 floats per material ID. palette[0] is unused (empty).
///          If null, all solid voxels get white (1,1,1,1).
/// filterMaterial: if >= 0, only that material ID is treated as solid; all
/// others read as empty. Avoids a separate JS-side masking pass.
///
/// borderX/Y/Z: per-axis halo size. The grid is expected to include `borderN`
/// extra voxels on each side of axis N (total dim = interior + 2*borderN).
/// Those halo voxels contribute to visibility queries (so faces facing a
/// solid halo voxel are suppressed), but no faces are emitted for them and
/// output positions are expressed in interior coordinates (the first interior
/// voxel sits at world-space origin). This lets chunked terrain pass in a
/// 1-voxel ring of neighbour data so shared chunk boundaries don't double up.
MeshData greedyMesh(const uint8_t* voxels, int gridX, int gridY, int gridZ,
                    float cellSize = 1.0f, const float* palette = nullptr,
                    int paletteCount = 0, int filterMaterial = -1,
                    int borderX = 0, int borderY = 0, int borderZ = 0);

} // namespace bromesh
