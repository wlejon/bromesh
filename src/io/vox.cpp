#include "bromesh/io/vox.h"

#include <cstdio>
#include <cstdint>
#include <cstring>

namespace bromesh {

VoxData loadVOX(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return {};

    // Read magic "VOX "
    char magic[4];
    if (std::fread(magic, 1, 4, f) != 4) { std::fclose(f); return {}; }
    if (std::memcmp(magic, "VOX ", 4) != 0) { std::fclose(f); return {}; }

    // Read version
    int32_t version = 0;
    if (std::fread(&version, sizeof(int32_t), 1, f) != 1) { std::fclose(f); return {}; }

    VoxData data;
    bool hasPalette = false;

    // Temporary storage for XYZI data (need SIZE first to allocate grid)
    struct VoxelEntry { uint8_t x, y, z, colorIndex; };
    std::vector<VoxelEntry> voxelEntries;

    // Read chunks
    // The first chunk is MAIN, which contains all other chunks as children.
    // We just read linearly through the file.
    auto readChunk = [&](auto& self) -> bool {
        char chunkId[4];
        if (std::fread(chunkId, 1, 4, f) != 4) return false;

        int32_t contentSize = 0, childrenSize = 0;
        if (std::fread(&contentSize, sizeof(int32_t), 1, f) != 1) return false;
        if (std::fread(&childrenSize, sizeof(int32_t), 1, f) != 1) return false;

        if (std::memcmp(chunkId, "MAIN", 4) == 0) {
            // MAIN chunk has no content, only children. Read children.
            long childEnd = std::ftell(f) + childrenSize;
            while (std::ftell(f) < childEnd) {
                if (!self(self)) break;
            }
        } else if (std::memcmp(chunkId, "SIZE", 4) == 0) {
            int32_t sx, sy, sz;
            if (std::fread(&sx, sizeof(int32_t), 1, f) != 1) return false;
            if (std::fread(&sy, sizeof(int32_t), 1, f) != 1) return false;
            if (std::fread(&sz, sizeof(int32_t), 1, f) != 1) return false;
            data.sizeX = sx;
            data.sizeY = sy;
            data.sizeZ = sz;
            // Skip remaining content if any
            int remaining = contentSize - 12;
            if (remaining > 0) std::fseek(f, remaining, SEEK_CUR);
        } else if (std::memcmp(chunkId, "XYZI", 4) == 0) {
            int32_t numVoxels = 0;
            if (std::fread(&numVoxels, sizeof(int32_t), 1, f) != 1) return false;
            voxelEntries.resize(numVoxels);
            for (int32_t i = 0; i < numVoxels; ++i) {
                uint8_t xyzi[4];
                if (std::fread(xyzi, 1, 4, f) != 4) return false;
                voxelEntries[i] = { xyzi[0], xyzi[1], xyzi[2], xyzi[3] };
            }
            // Skip remaining content if any
            int remaining = contentSize - 4 - numVoxels * 4;
            if (remaining > 0) std::fseek(f, remaining, SEEK_CUR);
        } else if (std::memcmp(chunkId, "RGBA", 4) == 0) {
            // 256 palette entries (but index 0 is unused in MagicaVoxel;
            // the file stores entries 1-255 then a padding entry)
            uint8_t rgba[256 * 4];
            size_t toRead = (contentSize < (int)sizeof(rgba)) ? contentSize : sizeof(rgba);
            if (std::fread(rgba, 1, toRead, f) != toRead) return false;
            // MagicaVoxel palette: indices 1-255 are stored as entries 0-254,
            // entry 255 is unused. Map file entry i -> palette index i+1.
            for (int i = 0; i < 255; ++i) {
                int palIdx = i + 1;
                data.palette[palIdx * 4 + 0] = rgba[i * 4 + 0] / 255.0f;
                data.palette[palIdx * 4 + 1] = rgba[i * 4 + 1] / 255.0f;
                data.palette[palIdx * 4 + 2] = rgba[i * 4 + 2] / 255.0f;
                data.palette[palIdx * 4 + 3] = rgba[i * 4 + 3] / 255.0f;
            }
            hasPalette = true;
            // Skip remaining content
            int remaining = contentSize - (int)toRead;
            if (remaining > 0) std::fseek(f, remaining, SEEK_CUR);
        } else {
            // Unknown chunk, skip content + children
            std::fseek(f, contentSize + childrenSize, SEEK_CUR);
        }
        return true;
    };

    readChunk(readChunk);

    std::fclose(f);

    // If no palette was found, generate a default one
    if (!hasPalette) {
        for (int i = 1; i < 256; ++i) {
            data.palette[i * 4 + 0] = 1.0f;
            data.palette[i * 4 + 1] = 1.0f;
            data.palette[i * 4 + 2] = 1.0f;
            data.palette[i * 4 + 3] = 1.0f;
        }
    }

    // Allocate and fill voxel grid
    if (data.sizeX > 0 && data.sizeY > 0 && data.sizeZ > 0) {
        size_t gridSize = (size_t)data.sizeX * data.sizeY * data.sizeZ;
        data.voxels.resize(gridSize, 0);

        for (auto& ve : voxelEntries) {
            if (ve.x < data.sizeX && ve.y < data.sizeY && ve.z < data.sizeZ) {
                size_t idx = (size_t)ve.z * data.sizeX * data.sizeY + (size_t)ve.y * data.sizeX + ve.x;
                data.voxels[idx] = ve.colorIndex;
            }
        }
    }

    return data;
}

} // namespace bromesh
