#include "test_framework.h"
#include <cmath>

static void fillSphereField(float* field, int N, float radius) {
    float c = (N - 1) * 0.5f;
    for (int z = 0; z < N; ++z)
        for (int y = 0; y < N; ++y)
            for (int x = 0; x < N; ++x) {
                float dx = x - c, dy = y - c, dz = z - c;
                field[z * N * N + y * N + x] = std::sqrt(dx*dx + dy*dy + dz*dz) - radius;
            }
}

TEST(reconstruct_from_sphere_pointcloud) {
    // Create a point cloud from a sphere
    auto sphere = bromesh::sphere(1.0f, 16, 12);
    bromesh::computeNormals(sphere);

    bromesh::ReconstructParams params;
    params.gridResolution = 32;

    auto result = bromesh::reconstructFromPointCloud(sphere, params);
    ASSERT(!result.empty(), "reconstruct: should produce mesh");
    ASSERT(result.triangleCount() > 10, "reconstruct: should have reasonable triangle count");
}

