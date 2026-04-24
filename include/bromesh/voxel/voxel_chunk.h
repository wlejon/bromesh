#pragma once

#include "bromesh/mesh_data.h"

#include <cstdint>
#include <vector>

namespace bromesh {

/// Fixed-height voxel chunk. Grid layout: voxels[(z * sizeY + y) * sizeX + x]
/// matching the greedyMesh convention.
class VoxelChunk {
public:
    /// Construct a chunk with given dimensions. All voxels init to 0 (air).
    VoxelChunk(int sizeX, int sizeY, int sizeZ, float cellSize = 1.0f);

    int sizeX() const { return sizeX_; }
    int sizeY() const { return sizeY_; }
    int sizeZ() const { return sizeZ_; }
    float cellSize() const { return cellSize_; }

    /// Get a voxel. Out-of-bounds returns 0.
    uint8_t getVoxel(int x, int y, int z) const;

    /// Set a voxel. Out-of-bounds is a no-op. Marks dirty.
    void setVoxel(int x, int y, int z, uint8_t material);

    /// Direct pointer to the voxel data for bulk fill (e.g. from noise).
    /// Caller must call markDirty() after modifying data directly.
    uint8_t* data() { return voxels_.data(); }
    const uint8_t* data() const { return voxels_.data(); }

    /// Fill all voxels with a value (0 to clear). Marks dirty.
    void fill(uint8_t value);

    bool isDirty() const { return dirty_; }
    void markDirty() { dirty_ = true; }
    void clearDirty() { dirty_ = false; }

    /// Build a greedy mesh with per-vertex colors from palette.
    /// All materials in one pass (filterMaterial = -1).
    /// palette: RGBA floats, 4 per material ID. Index 0 unused (air).
    /// borderX/Y/Z: halo voxels on each side that participate in visibility
    /// checks but don't produce faces. Output positions are expressed in
    /// interior coordinates (first interior voxel at local origin). See
    /// `greedyMesh` for the full semantics.
    /// Returns empty MeshData if chunk has no solid voxels.
    MeshData buildMesh(const float* palette = nullptr, int paletteCount = 0,
                       int borderX = 0, int borderY = 0, int borderZ = 0) const;

private:
    bool inBounds(int x, int y, int z) const {
        return x >= 0 && x < sizeX_ &&
               y >= 0 && y < sizeY_ &&
               z >= 0 && z < sizeZ_;
    }

    int index(int x, int y, int z) const {
        return (z * sizeY_ + y) * sizeX_ + x;
    }

    int sizeX_, sizeY_, sizeZ_;
    float cellSize_;
    std::vector<uint8_t> voxels_;
    bool dirty_ = true;
};

} // namespace bromesh
