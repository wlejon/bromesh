#include "bromesh/analysis/sdf.h"
#include "bromesh/analysis/bvh.h"
#include "bromesh/analysis/bbox.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace bromesh {

namespace {

int countRayCrossings(const MeshBVH& bvh, const MeshData& mesh,
                      const float origin[3], const float dir[3]) {
    int crossings = 0;
    float curOrigin[3] = { origin[0], origin[1], origin[2] };
    for (int iter = 0; iter < 64; ++iter) {
        RayHit rh = bvh.raycast(mesh, curOrigin, dir);
        if (!rh.hit) break;
        crossings++;
        constexpr float kEps = 1e-4f;
        curOrigin[0] = rh.position[0] + dir[0] * kEps;
        curOrigin[1] = rh.position[1] + dir[1] * kEps;
        curOrigin[2] = rh.position[2] + dir[2] * kEps;
    }
    return crossings;
}

} // namespace

SDFVolume generateSDF(const MeshData& mesh, int dimX, int dimY, int dimZ,
                      const SDFOptions& opts) {
    SDFVolume vol;
    if (dimX <= 0 || dimY <= 0 || dimZ <= 0) {
        return vol;
    }
    vol.dimX = dimX;
    vol.dimY = dimY;
    vol.dimZ = dimZ;

    if (mesh.empty() || mesh.indices.empty()) {
        vol.field.assign(static_cast<size_t>(dimX) * dimY * dimZ, 1e9f);
        return vol;
    }

    bool boundsValid = !bromath::aisEmpty(opts.bounds) &&
                       (opts.bounds.max.x > opts.bounds.min.x) &&
                       (opts.bounds.max.y > opts.bounds.min.y) &&
                       (opts.bounds.max.z > opts.bounds.min.z);
    bromath::AABB3 bounds;
    if (boundsValid) {
        bounds = opts.bounds;
    } else {
        bromath::AABB3 meshBounds = computeBBox(mesh);
        bromath::Vec3 ext = bromath::aextent(meshBounds);
        float maxExtent = std::max({ext.x, ext.y, ext.z, 1.0f});
        float pad = opts.padding * maxExtent;
        bounds.min = meshBounds.min - bromath::Vec3{pad, pad, pad};
        bounds.max = meshBounds.max + bromath::Vec3{pad, pad, pad};
    }
    vol.bounds = bounds;

    float cellSize[3] = {
        dimX > 1 ? (bounds.max.x - bounds.min.x) / static_cast<float>(dimX - 1) : 1.0f,
        dimY > 1 ? (bounds.max.y - bounds.min.y) / static_cast<float>(dimY - 1) : 1.0f,
        dimZ > 1 ? (bounds.max.z - bounds.min.z) / static_cast<float>(dimZ - 1) : 1.0f
    };
    vol.cellSize[0] = cellSize[0];
    vol.cellSize[1] = cellSize[1];
    vol.cellSize[2] = cellSize[2];

    MeshBVH bvh = MeshBVH::build(mesh);
    if (bvh.empty()) {
        vol.field.assign(static_cast<size_t>(dimX) * dimY * dimZ, 1e9f);
        return vol;
    }

    vol.field.resize(static_cast<size_t>(dimX) * dimY * dimZ, 0.0f);

    float rayDir[3] = { 0.99999f, 0.001337f, 0.002718f };
    float rayDirLen = std::sqrt(rayDir[0]*rayDir[0] + rayDir[1]*rayDir[1] + rayDir[2]*rayDir[2]);
    rayDir[0] /= rayDirLen; rayDir[1] /= rayDirLen; rayDir[2] /= rayDirLen;

    float altDir[3] = { 0.0017f, 0.99999f, 0.0023f };
    float altDirLen = std::sqrt(altDir[0]*altDir[0] + altDir[1]*altDir[1] + altDir[2]*altDir[2]);
    altDir[0] /= altDirLen; altDir[1] /= altDirLen; altDir[2] /= altDirLen;

    #pragma omp parallel for
    for (int z = 0; z < dimZ; ++z) {
        float pz = bounds.min.z + z * cellSize[2];
        for (int y = 0; y < dimY; ++y) {
            float py = bounds.min.y + y * cellSize[1];
            for (int x = 0; x < dimX; ++x) {
                float px = bounds.min.x + x * cellSize[0];
                float p[3] = { px, py, pz };

                RayHit hit = bvh.closestPoint(mesh, p);
                if (!hit.hit) {
                    vol.field[(z * dimY + y) * dimX + x] = 1e9f;
                    continue;
                }

                float diff[3] = {
                    px - hit.position[0],
                    py - hit.position[1],
                    pz - hit.position[2]
                };
                float dotVal = diff[0] * hit.normal[0] +
                               diff[1] * hit.normal[1] +
                               diff[2] * hit.normal[2];

                bool dotInside = (dotVal < 0.0f);
                int crossings = countRayCrossings(bvh, mesh, p, rayDir);
                bool rayInside = (crossings % 2) != 0;
                bool inside = rayInside;

                if (rayInside != dotInside) {
                    int altCrossings = countRayCrossings(bvh, mesh, p, altDir);
                    bool altInside = (altCrossings % 2) != 0;
                    if (altInside == rayInside) {
                        inside = rayInside;
                    } else if (altInside == dotInside) {
                        inside = dotInside;
                    } else {
                        inside = rayInside;
                    }
                }

                float sign = inside ? -1.0f : 1.0f;
                vol.field[(z * dimY + y) * dimX + x] = sign * hit.distance;
            }
        }
    }

    return vol;
}

} // namespace bromesh
