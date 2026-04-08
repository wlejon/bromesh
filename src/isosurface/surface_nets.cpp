#include "bromesh/isosurface/surface_nets.h"

namespace bromesh {

MeshData surfaceNets(const float* field, int gridX, int gridY, int gridZ,
                     float isoLevel, float cellSize) {
    (void)field; (void)gridX; (void)gridY; (void)gridZ;
    (void)isoLevel; (void)cellSize;
    return {};
}

} // namespace bromesh
