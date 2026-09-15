#include "test_framework.h"
#include "bromesh/analysis/sdf.h"
#include "bromesh/primitives/primitives.h"
#include "bromesh/isosurface/marching_cubes.h"
#include "bromesh/isosurface/dual_contouring.h"
#include "bromesh/analysis/bbox.h"

#include <cmath>
#include <algorithm>

TEST(sdf_sphere_analytic_agreement) {
    const float radius = 1.0f;
    auto sph = bromesh::sphere(radius, 64, 32);
    const int dim = 32;

    bromesh::SDFOptions opts;
    opts.padding = 0.2f;
    auto sdf = bromesh::generateSDF(sph, dim, dim, dim, opts);

    ASSERT(sdf.dimX == dim && sdf.dimY == dim && sdf.dimZ == dim, "grid dims match");
    ASSERT(sdf.field.size() == static_cast<size_t>(dim * dim * dim), "field size matches");

    float voxelSize = sdf.cellSize[0];
    float maxDiff = 0.0f;

    for (int z = 0; z < dim; ++z) {
        float pz = sdf.bounds.min.z + z * sdf.cellSize[2];
        for (int y = 0; y < dim; ++y) {
            float py = sdf.bounds.min.y + y * sdf.cellSize[1];
            for (int x = 0; x < dim; ++x) {
                float px = sdf.bounds.min.x + x * sdf.cellSize[0];
                float pLen = std::sqrt(px * px + py * py + pz * pz);
                float analyticSDF = pLen - radius;
                float computedSDF = sdf.value(x, y, z);
                float diff = std::fabs(computedSDF - analyticSDF);
                if (diff > maxDiff) {
                    maxDiff = diff;
                }
            }
        }
    }

    // Verify maximum deviation is small (< 1-2 voxels)
    ASSERT(maxDiff < 2.0f * voxelSize, "max deviation is < 2 voxels");
}

TEST(sdf_isosurface_round_trip) {
    bromesh::SDFOptions opts;
    opts.padding = 0.15f;

    // 1. Box test with marching cubes
    {
        auto b = bromesh::box(1.0f, 1.0f, 1.0f);
        auto origBBox = bromesh::computeBBox(b);

        auto sdfBox = bromesh::generateSDF(b, 32, 32, 32, opts);
        auto mcBox = bromesh::marchingCubes(sdfBox.field.data(),
                                            sdfBox.dimX, sdfBox.dimY, sdfBox.dimZ,
                                            0.0f, 1.0f);

        ASSERT(!mcBox.empty(), "extracted box non-empty");
        ASSERT(mcBox.triangleCount() > 0, "extracted box has triangles");

        for (size_t i = 0; i < mcBox.vertexCount(); ++i) {
            mcBox.positions[i * 3 + 0] = sdfBox.bounds.min.x + mcBox.positions[i * 3 + 0] * sdfBox.cellSize[0];
            mcBox.positions[i * 3 + 1] = sdfBox.bounds.min.y + mcBox.positions[i * 3 + 1] * sdfBox.cellSize[1];
            mcBox.positions[i * 3 + 2] = sdfBox.bounds.min.z + mcBox.positions[i * 3 + 2] * sdfBox.cellSize[2];
        }

        float volBox = bromesh::computeVolume(mcBox);
        ASSERT(volBox > 0.0f, "extracted box has positive volume");
        ASSERT(std::fabs(volBox - 8.0f) < 2.0f, "extracted box volume near original 8.0");

        auto mcBBox = bromesh::computeBBox(mcBox);
        float voxX = sdfBox.cellSize[0];
        float voxY = sdfBox.cellSize[1];
        float voxZ = sdfBox.cellSize[2];
        ASSERT(std::fabs(mcBBox.min.x - origBBox.min.x) < 2.0f * voxX, "box min.x matches original bounds");
        ASSERT(std::fabs(mcBBox.max.x - origBBox.max.x) < 2.0f * voxX, "box max.x matches original bounds");
        ASSERT(std::fabs(mcBBox.min.y - origBBox.min.y) < 2.0f * voxY, "box min.y matches original bounds");
        ASSERT(std::fabs(mcBBox.max.y - origBBox.max.y) < 2.0f * voxY, "box max.y matches original bounds");
        ASSERT(std::fabs(mcBBox.min.z - origBBox.min.z) < 2.0f * voxZ, "box min.z matches original bounds");
        ASSERT(std::fabs(mcBBox.max.z - origBBox.max.z) < 2.0f * voxZ, "box max.z matches original bounds");
    }

    // 2. Torus test with dual contouring / marching cubes
    {
        auto t = bromesh::torus(1.0f, 0.35f, 32, 16);
        auto origBBox = bromesh::computeBBox(t);

        auto sdfTorus = bromesh::generateSDF(t, 32, 32, 32, opts);
        auto mcTorus = bromesh::marchingCubes(sdfTorus.field.data(),
                                              sdfTorus.dimX, sdfTorus.dimY, sdfTorus.dimZ,
                                              0.0f, 1.0f);

        ASSERT(!mcTorus.empty(), "extracted torus non-empty");
        ASSERT(mcTorus.triangleCount() > 0, "extracted torus has triangles");

        for (size_t i = 0; i < mcTorus.vertexCount(); ++i) {
            mcTorus.positions[i * 3 + 0] = sdfTorus.bounds.min.x + mcTorus.positions[i * 3 + 0] * sdfTorus.cellSize[0];
            mcTorus.positions[i * 3 + 1] = sdfTorus.bounds.min.y + mcTorus.positions[i * 3 + 1] * sdfTorus.cellSize[1];
            mcTorus.positions[i * 3 + 2] = sdfTorus.bounds.min.z + mcTorus.positions[i * 3 + 2] * sdfTorus.cellSize[2];
        }

        float volTorus = bromesh::computeVolume(mcTorus);
        ASSERT(volTorus > 0.0f, "extracted torus has positive volume");

        auto mcBBox = bromesh::computeBBox(mcTorus);
        float voxX = sdfTorus.cellSize[0];
        float voxY = sdfTorus.cellSize[1];
        float voxZ = sdfTorus.cellSize[2];
        ASSERT(std::fabs(mcBBox.min.x - origBBox.min.x) < 2.0f * voxX, "torus min.x matches original bounds");
        ASSERT(std::fabs(mcBBox.max.x - origBBox.max.x) < 2.0f * voxX, "torus max.x matches original bounds");
        ASSERT(std::fabs(mcBBox.min.y - origBBox.min.y) < 2.0f * voxY, "torus min.y matches original bounds");
        ASSERT(std::fabs(mcBBox.max.y - origBBox.max.y) < 2.0f * voxY, "torus max.y matches original bounds");
        ASSERT(std::fabs(mcBBox.min.z - origBBox.min.z) < 2.0f * voxZ, "torus min.z matches original bounds");
        ASSERT(std::fabs(mcBBox.max.z - origBBox.max.z) < 2.0f * voxZ, "torus max.z matches original bounds");
    }
}
