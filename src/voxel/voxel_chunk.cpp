#include "bromesh/voxel/voxel_chunk.h"
#include "bromesh/voxel/greedy_mesh.h"

#include <algorithm>
#include <cstring>

namespace bromesh {

VoxelChunk::VoxelChunk(int sizeX, int sizeY, int sizeZ, float cellSize)
    : sizeX_(sizeX), sizeY_(sizeY), sizeZ_(sizeZ), cellSize_(cellSize),
      voxels_(static_cast<size_t>(sizeX) * sizeY * sizeZ, 0) {}

uint8_t VoxelChunk::getVoxel(int x, int y, int z) const {
    if (!inBounds(x, y, z)) return 0;
    return voxels_[index(x, y, z)];
}

void VoxelChunk::setVoxel(int x, int y, int z, uint8_t material) {
    if (!inBounds(x, y, z)) return;
    voxels_[index(x, y, z)] = material;
    dirty_ = true;
}

void VoxelChunk::fill(uint8_t value) {
    std::memset(voxels_.data(), value, voxels_.size());
    dirty_ = true;
}

MeshData VoxelChunk::buildMesh(const float* palette, int paletteCount,
                               int borderX, int borderY, int borderZ) const {
    return greedyMesh(voxels_.data(), sizeX_, sizeY_, sizeZ_, cellSize_,
                      palette, paletteCount, /*filterMaterial=*/-1,
                      borderX, borderY, borderZ);
}

} // namespace bromesh
