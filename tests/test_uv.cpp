#include "test_framework.h"
#include <cmath>

TEST(project_uvs_box) {
    auto b = bromesh::box(1, 1, 1);
    b.uvs.clear(); // remove existing UVs
    ASSERT(!b.hasUVs(), "UVs should be cleared");

    bromesh::projectUVs(b, bromesh::ProjectionType::Box, 1.0f);
    ASSERT(b.hasUVs(), "box projection should produce UVs");
    ASSERT(b.uvs.size() == b.vertexCount() * 2, "UV count matches vertex count");

    // Test planar XY too
    auto b2 = bromesh::box(1, 1, 1);
    b2.uvs.clear();
    bromesh::projectUVs(b2, bromesh::ProjectionType::PlanarXY, 1.0f);
    ASSERT(b2.hasUVs(), "planar XY projection should produce UVs");
}


#if BROMESH_HAS_XATLAS
TEST(unwrap_sphere) {
    auto mesh = bromesh::sphere(2.0f, 16, 12);
    bromesh::computeNormals(mesh);
    size_t origTriCount = mesh.triangleCount();

    auto result = bromesh::unwrapUVs(mesh);
    ASSERT(result.success, "unwrap_sphere: should succeed");
    ASSERT(result.chartCount > 0, "unwrap_sphere: should have charts");
    ASSERT(result.atlasWidth > 0, "unwrap_sphere: atlas width > 0");
    ASSERT(result.atlasHeight > 0, "unwrap_sphere: atlas height > 0");
    ASSERT(mesh.hasUVs(), "unwrap_sphere: should have UVs after unwrap");
    ASSERT(mesh.triangleCount() == origTriCount, "unwrap_sphere: triangle count preserved");

    // UVs should be in [0,1] range
    bool uvsInRange = true;
    for (size_t i = 0; i < mesh.uvs.size(); ++i) {
        if (mesh.uvs[i] < -0.01f || mesh.uvs[i] > 1.01f) {
            uvsInRange = false; break;
        }
    }
    ASSERT(uvsInRange, "unwrap_sphere: UVs should be in [0,1] range");
}

TEST(unwrap_box) {
    auto mesh = bromesh::box(1, 1, 1);
    bromesh::computeNormals(mesh);

    auto result = bromesh::unwrapUVs(mesh);
    ASSERT(result.success, "unwrap_box: should succeed");
    ASSERT(mesh.hasUVs(), "unwrap_box: should have UVs");
    ASSERT(mesh.hasNormals(), "unwrap_box: should preserve normals");
}

TEST(unwrap_torus) {
    auto mesh = bromesh::torus(2.0f, 0.5f, 24, 12);
    bromesh::computeNormals(mesh);

    auto result = bromesh::unwrapUVs(mesh);
    ASSERT(result.success, "unwrap_torus: should succeed");
    ASSERT(result.chartCount > 1, "unwrap_torus: torus should have multiple charts");
    ASSERT(mesh.hasUVs(), "unwrap_torus: should have UVs");
}

TEST(unwrap_custom_params) {
    auto mesh = bromesh::cylinder(1.0f, 2.0f, 16);
    bromesh::computeNormals(mesh);

    bromesh::UnwrapParams cp;
    cp.maxStretch = 0.1f;
    bromesh::PackParams pp;
    pp.padding = 2;

    auto result = bromesh::unwrapUVs(mesh, cp, pp);
    ASSERT(result.success, "unwrap_custom: should succeed");
    ASSERT(mesh.hasUVs(), "unwrap_custom: should have UVs");
}

#endif // BROMESH_HAS_XATLAS

TEST(uv_metrics_box_projection) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Box, 1.0f);

    auto metrics = bromesh::measureUVQuality(mesh);
    ASSERT(metrics.triangleCount == mesh.triangleCount(),
           "uv_metrics: triangle count matches");
    ASSERT(metrics.avgStretch > 0.0f, "uv_metrics: avg stretch > 0");
    ASSERT(metrics.uvSpaceUsage > 0.0f, "uv_metrics: UV usage > 0");
}

TEST(uv_distortion_per_triangle) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Box, 1.0f);

    auto distortions = bromesh::computeUVDistortion(mesh);
    ASSERT(distortions.size() == mesh.triangleCount(),
           "uv_distortion: one entry per triangle");

    // All stretch values should be positive
    bool allPositive = true;
    for (auto& d : distortions) {
        if (d.stretch <= 0.0f) { allPositive = false; break; }
    }
    ASSERT(allPositive, "uv_distortion: all stretch values positive");
}

TEST(uv_metrics_no_uvs) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    mesh.uvs.clear(); // remove UVs

    auto metrics = bromesh::measureUVQuality(mesh);
    ASSERT(metrics.triangleCount == 0, "uv_metrics_no_uvs: no triangles when no UVs");

    auto distortions = bromesh::computeUVDistortion(mesh);
    ASSERT(distortions.empty(), "uv_distortion_no_uvs: empty when no UVs");
}

TEST(uv_metrics_planar_plane) {
    // A flat plane with planar XZ projection should have low angle distortion
    auto mesh = bromesh::plane(2.0f, 2.0f, 4, 4);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::PlanarXZ, 1.0f);

    auto metrics = bromesh::measureUVQuality(mesh);
    ASSERT(metrics.avgAngleDistortion < 0.1f,
           "uv_planar: flat plane should have low angle distortion");
}

