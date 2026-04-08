#include "bromesh/voxel/greedy_mesh.h"

namespace bromesh {

MeshData greedyMesh(const uint8_t* voxels, int gridX, int gridY, int gridZ,
                    float cellSize, const float* palette, int paletteCount) {
    (void)voxels; (void)gridX; (void)gridY; (void)gridZ;
    (void)cellSize; (void)palette; (void)paletteCount;
    return {};
}

} // namespace bromesh
