#pragma once
#include "bromesh/mesh_data.h"
#include <bromath/aabb.h>
#include <vector>

namespace bromesh {

struct SDFVolume {
    std::vector<float> field; // row-major: (z * dimY + y) * dimX + x
    int dimX = 0;
    int dimY = 0;
    int dimZ = 0;
    bromath::AABB3 bounds;
    float cellSize[3] = {1.0f, 1.0f, 1.0f};

    float value(int x, int y, int z) const {
        if (x < 0 || x >= dimX || y < 0 || y >= dimY || z < 0 || z >= dimZ)
            return 1e9f;
        return field[(z * dimY + y) * dimX + x];
    }
};

struct SDFOptions {
    /// Optional explicit bounds. If empty/inverted, computed from mesh bounds with padding.
    bromath::AABB3 bounds;
    /// Extra padding fraction around mesh bounding box (default 0.1 = 10%).
    float padding = 0.1f;
};

/// Compute a signed distance field (SDF) on a uniform 3D grid from a closed surface mesh.
/// Sign convention follows bromesh isosurface: < 0 inside, >= 0 outside.
/// Accelerated via MeshBVH::closestPoint and ray-parity winding.
SDFVolume generateSDF(const MeshData& mesh, int dimX, int dimY, int dimZ,
                      const SDFOptions& opts = {});

} // namespace bromesh
