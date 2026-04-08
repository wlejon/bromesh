#include "bromesh/isosurface/transvoxel.h"

namespace bromesh {

MeshData transvoxel(const float* field, int gridSize, int lod,
                    const int neighborLods[6],
                    float isoLevel, float cellSize) {
    (void)field; (void)gridSize; (void)lod; (void)neighborLods;
    (void)isoLevel; (void)cellSize;
    return {};
}

} // namespace bromesh
