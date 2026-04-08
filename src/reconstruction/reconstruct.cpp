#include "bromesh/reconstruction/reconstruct.h"
#include "bromesh/isosurface/marching_cubes.h"
#include "bromesh/analysis/bbox.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

namespace bromesh {

// Implicit surface reconstruction via signed-distance estimation:
// For each voxel, estimate the signed distance to the surface by
// accumulating weighted projections onto point normals.

MeshData reconstructFromPointCloud(
    const float* positions, const float* normals, size_t pointCount,
    const ReconstructParams& params) {

    if (!positions || !normals || pointCount == 0) return {};

    // Compute bounding box with padding
    float bmin[3] = { positions[0], positions[1], positions[2] };
    float bmax[3] = { positions[0], positions[1], positions[2] };

    for (size_t i = 1; i < pointCount; ++i) {
        for (int c = 0; c < 3; ++c) {
            bmin[c] = std::min(bmin[c], positions[i * 3 + c]);
            bmax[c] = std::max(bmax[c], positions[i * 3 + c]);
        }
    }

    float extent[3] = { bmax[0]-bmin[0], bmax[1]-bmin[1], bmax[2]-bmin[2] };
    float maxExtent = std::max({extent[0], extent[1], extent[2]});
    if (maxExtent < 1e-8f) return {};

    // Add 10% padding
    float pad = maxExtent * 0.1f;
    for (int c = 0; c < 3; ++c) {
        bmin[c] -= pad;
        bmax[c] += pad;
        extent[c] = bmax[c] - bmin[c];
    }

    int res = std::max(8, params.gridResolution);
    float cellSize = maxExtent / (res - 1);

    int gridX = static_cast<int>(std::ceil(extent[0] / cellSize)) + 1;
    int gridY = static_cast<int>(std::ceil(extent[1] / cellSize)) + 1;
    int gridZ = static_cast<int>(std::ceil(extent[2] / cellSize)) + 1;

    // Clamp to avoid excessive memory
    gridX = std::min(gridX, 256);
    gridY = std::min(gridY, 256);
    gridZ = std::min(gridZ, 256);

    // Support radius: auto-estimate from average nearest-neighbor distance
    float supportRadius = params.supportRadius;
    if (supportRadius <= 0.0f) {
        // Heuristic: assume roughly uniform distribution
        float volume = extent[0] * extent[1] * extent[2];
        float avgSpacing = std::cbrt(volume / pointCount);
        supportRadius = avgSpacing * 3.0f;
    }
    float radiusSq = supportRadius * supportRadius;

    // Build spatial hash for fast neighbor queries
    float invCell = 1.0f / supportRadius;
    struct HashEntry { size_t pointIdx; };

    // Simple grid-based spatial lookup
    int hashDimX = std::max(1, static_cast<int>(extent[0] * invCell) + 1);
    int hashDimY = std::max(1, static_cast<int>(extent[1] * invCell) + 1);
    int hashDimZ = std::max(1, static_cast<int>(extent[2] * invCell) + 1);

    // Clamp hash dimensions
    hashDimX = std::min(hashDimX, 128);
    hashDimY = std::min(hashDimY, 128);
    hashDimZ = std::min(hashDimZ, 128);

    size_t hashSize = static_cast<size_t>(hashDimX) * hashDimY * hashDimZ;
    std::vector<std::vector<size_t>> hashGrid(hashSize);

    for (size_t i = 0; i < pointCount; ++i) {
        int hx = std::min(hashDimX - 1, std::max(0, static_cast<int>((positions[i*3+0] - bmin[0]) * invCell)));
        int hy = std::min(hashDimY - 1, std::max(0, static_cast<int>((positions[i*3+1] - bmin[1]) * invCell)));
        int hz = std::min(hashDimZ - 1, std::max(0, static_cast<int>((positions[i*3+2] - bmin[2]) * invCell)));
        size_t idx = static_cast<size_t>(hx) + static_cast<size_t>(hy) * hashDimX
            + static_cast<size_t>(hz) * hashDimX * hashDimY;
        hashGrid[idx].push_back(i);
    }

    // Evaluate implicit function on voxel grid
    size_t fieldSize = static_cast<size_t>(gridX) * gridY * gridZ;
    std::vector<float> field(fieldSize, 1.0f); // positive = outside

    for (int iz = 0; iz < gridZ; ++iz) {
        for (int iy = 0; iy < gridY; ++iy) {
            for (int ix = 0; ix < gridX; ++ix) {
                float px = bmin[0] + ix * cellSize;
                float py = bmin[1] + iy * cellSize;
                float pz = bmin[2] + iz * cellSize;

                // Find hash cell
                int hx = std::min(hashDimX - 1, std::max(0, static_cast<int>((px - bmin[0]) * invCell)));
                int hy = std::min(hashDimY - 1, std::max(0, static_cast<int>((py - bmin[1]) * invCell)));
                int hz = std::min(hashDimZ - 1, std::max(0, static_cast<int>((pz - bmin[2]) * invCell)));

                float weightedDist = 0.0f;
                float totalWeight = 0.0f;

                // Search 3x3x3 neighborhood
                for (int dz = -1; dz <= 1; ++dz) {
                    int nz = hz + dz;
                    if (nz < 0 || nz >= hashDimZ) continue;
                    for (int dy = -1; dy <= 1; ++dy) {
                        int ny = hy + dy;
                        if (ny < 0 || ny >= hashDimY) continue;
                        for (int dx = -1; dx <= 1; ++dx) {
                            int nx = hx + dx;
                            if (nx < 0 || nx >= hashDimX) continue;

                            size_t hidx = static_cast<size_t>(nx)
                                + static_cast<size_t>(ny) * hashDimX
                                + static_cast<size_t>(nz) * hashDimX * hashDimY;

                            for (size_t pi : hashGrid[hidx]) {
                                float dpx = px - positions[pi*3+0];
                                float dpy = py - positions[pi*3+1];
                                float dpz = pz - positions[pi*3+2];
                                float distSq = dpx*dpx + dpy*dpy + dpz*dpz;

                                if (distSq < radiusSq) {
                                    // Wendland C2 kernel
                                    float r = std::sqrt(distSq) / supportRadius;
                                    float h = 1.0f - r;
                                    float w = h * h * h * h * (4.0f * r + 1.0f);

                                    // Signed distance: projection onto normal
                                    float sd = dpx * normals[pi*3+0]
                                             + dpy * normals[pi*3+1]
                                             + dpz * normals[pi*3+2];

                                    weightedDist += w * sd;
                                    totalWeight += w;
                                }
                            }
                        }
                    }
                }

                size_t vi = static_cast<size_t>(ix)
                    + static_cast<size_t>(iy) * gridX
                    + static_cast<size_t>(iz) * gridX * gridY;

                if (totalWeight > 1e-8f) {
                    field[vi] = weightedDist / totalWeight;
                } else {
                    field[vi] = supportRadius; // Far from surface
                }
            }
        }
    }

    // Extract surface with marching cubes
    MeshData result = marchingCubes(field.data(), gridX, gridY, gridZ, 0.0f, cellSize);

    // Offset positions to world space
    for (size_t v = 0; v < result.vertexCount(); ++v) {
        result.positions[v * 3 + 0] += bmin[0];
        result.positions[v * 3 + 1] += bmin[1];
        result.positions[v * 3 + 2] += bmin[2];
    }

    return result;
}

MeshData reconstructFromPointCloud(
    const MeshData& pointCloud,
    const ReconstructParams& params) {

    if (pointCloud.positions.empty() || !pointCloud.hasNormals()) return {};

    return reconstructFromPointCloud(
        pointCloud.positions.data(),
        pointCloud.normals.data(),
        pointCloud.vertexCount(),
        params);
}

} // namespace bromesh
