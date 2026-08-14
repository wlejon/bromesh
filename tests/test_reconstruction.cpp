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
    // Create a point cloud from a sphere of radius R=1.0
    auto sphere = bromesh::sphere(1.0f, 32, 24);
    bromesh::computeNormals(sphere);

    bromesh::ReconstructParams params;
    params.gridResolution = 32;
    params.supportRadius = 0.30f;

    auto rawResult = bromesh::reconstructFromPointCloud(sphere, params);
    ASSERT(!rawResult.empty(), "reconstruct: sphere should produce non-empty mesh");

    // Extract outer reconstructed surface component
    auto comps = bromesh::splitConnectedComponents(rawResult);
    ASSERT(!comps.empty(), "reconstruct: should have at least one component");
    // Select largest component by triangle count (outer surface)
    size_t bestIdx = 0;
    for (size_t i = 1; i < comps.size(); ++i) {
        if (comps[i].triangleCount() > comps[bestIdx].triangleCount()) {
            bestIdx = i;
        }
    }
    const auto& result = comps[bestIdx];
    ASSERT(!result.empty(), "reconstruct sphere: outer component non-empty");
    ASSERT(result.triangleCount() > 10, "reconstruct sphere: should have reasonable triangle count");

    // Assert every vertex p satisfies |||p|| - 1.0| < 0.15
    bool allVertsNearSphere = true;
    for (size_t v = 0; v < result.vertexCount(); ++v) {
        float px = result.positions[v * 3 + 0];
        float py = result.positions[v * 3 + 1];
        float pz = result.positions[v * 3 + 2];
        float dist = std::sqrt(px * px + py * py + pz * pz);
        if (std::fabs(dist - 1.0f) >= 0.15f) {
            allVertsNearSphere = false;
            break;
        }
    }
    ASSERT(allVertsNearSphere, "reconstruct sphere: all vertices satisfy ||p|| - 1.0 < 0.15");

    // Assert computeVolume(result) is close to analytic sphere volume 4/3*pi ~= 4.18879 (within 10%)
    float analyticSphereVol = (4.0f / 3.0f) * 3.14159265f;
    float resultVol = bromesh::computeVolume(result);
    ASSERT(std::fabs(resultVol - analyticSphereVol) < 0.10f * analyticSphereVol,
           "reconstruct sphere: volume matches 4/3*pi within 10%");

    // Assert bounding box extent is within [0.85, 1.15] * 2.0 = [1.70, 2.30]
    auto bbox = bromesh::computeBBox(result);
    float extX = bbox.max.x - bbox.min.x;
    float extY = bbox.max.y - bbox.min.y;
    float extZ = bbox.max.z - bbox.min.z;
    ASSERT(extX >= 1.70f && extX <= 2.30f, "reconstruct sphere: X extent within [0.85, 1.15]*2.0");
    ASSERT(extY >= 1.70f && extY <= 2.30f, "reconstruct sphere: Y extent within [0.85, 1.15]*2.0");
    ASSERT(extZ >= 1.70f && extZ <= 2.30f, "reconstruct sphere: Z extent within [0.85, 1.15]*2.0");
}

TEST(reconstruct_from_box_pointcloud) {
    // Sample points on a box of size 2x2x2 (half-extent 1.0, volume 8.0, bounds [-1, 1]^3)
    auto boxMesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::computeNormals(boxMesh);

    // Sample dense surface points (or generate a grid of points on each face)
    auto pointCloud = bromesh::sampleSurface(boxMesh, 2000, 12345);

    bromesh::ReconstructParams params;
    params.gridResolution = 32;

    auto result = bromesh::reconstructFromPointCloud(pointCloud, params);
    ASSERT(!result.empty(), "reconstruct box: should produce non-empty mesh");
    ASSERT(result.triangleCount() > 10, "reconstruct box: should have triangles");

    // Volume matches 8.0 +- 1.0
    float resultVol = bromesh::computeVolume(result);
    ASSERT(std::fabs(resultVol - 8.0f) <= 1.0f,
           "reconstruct box: volume matches 8.0 +- 1.0");

    // Bounding box matches [-1, 1]^3 +- 0.15
    auto bbox = bromesh::computeBBox(result);
    ASSERT(std::fabs(bbox.min.x - (-1.0f)) <= 0.15f && std::fabs(bbox.max.x - 1.0f) <= 0.15f,
           "reconstruct box: X bounds match [-1, 1] +- 0.15");
    ASSERT(std::fabs(bbox.min.y - (-1.0f)) <= 0.15f && std::fabs(bbox.max.y - 1.0f) <= 0.15f,
           "reconstruct box: Y bounds match [-1, 1] +- 0.15");
    ASSERT(std::fabs(bbox.min.z - (-1.0f)) <= 0.15f && std::fabs(bbox.max.z - 1.0f) <= 0.15f,
           "reconstruct box: Z bounds match [-1, 1] +- 0.15");
}


