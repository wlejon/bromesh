#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace bromesh {

/// MagicaVoxel .vox file data.
struct VoxData {
    int sizeX = 0, sizeY = 0, sizeZ = 0;
    std::vector<uint8_t> voxels; // 3D grid, 0 = empty, nonzero = palette index
    float palette[256 * 4] = {}; // RGBA float, 256 entries
};

/// Load a MagicaVoxel .vox file. Returns empty data on failure.
VoxData loadVOX(const std::string& path);

} // namespace bromesh
