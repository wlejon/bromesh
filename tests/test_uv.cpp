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

TEST(project_uvs_planar_xy_analytic) {
    // Quad with vertices at (-1, -1, 0), (1, -1, 0), (1, 1, 0), (-1, 1, 0)
    bromesh::MeshData mesh;
    mesh.positions = {
        -1.0f, -1.0f, 0.0f,
         1.0f, -1.0f, 0.0f,
         1.0f,  1.0f, 0.0f,
        -1.0f,  1.0f, 0.0f
    };
    mesh.indices = { 0, 1, 2, 0, 2, 3 };

    bromesh::projectUVs(mesh, bromesh::ProjectionType::PlanarXY, 1.0f);
    ASSERT(mesh.hasUVs(), "project_uvs_planar_xy_analytic: has UVs");
    ASSERT(std::fabs(mesh.uvs[0] - (-1.0f)) < 1e-4f && std::fabs(mesh.uvs[1] - (-1.0f)) < 1e-4f,
           "planar XY: (-1,-1) -> (-1,-1)");
    ASSERT(std::fabs(mesh.uvs[2] - (1.0f)) < 1e-4f && std::fabs(mesh.uvs[3] - (-1.0f)) < 1e-4f,
           "planar XY: (1,-1) -> (1,-1)");
    ASSERT(std::fabs(mesh.uvs[4] - (1.0f)) < 1e-4f && std::fabs(mesh.uvs[5] - (1.0f)) < 1e-4f,
           "planar XY: (1,1) -> (1,1)");
    ASSERT(std::fabs(mesh.uvs[6] - (-1.0f)) < 1e-4f && std::fabs(mesh.uvs[7] - (1.0f)) < 1e-4f,
           "planar XY: (-1,1) -> (-1,1)");

    // Unit quad in [0, 1] interval
    bromesh::MeshData unitQuad;
    unitQuad.positions = {
        0.0f, 0.0f, 0.0f,
        1.0f, 0.0f, 0.0f,
        1.0f, 1.0f, 0.0f,
        0.0f, 1.0f, 0.0f
    };
    unitQuad.indices = { 0, 1, 2, 0, 2, 3 };
    bromesh::projectUVs(unitQuad, bromesh::ProjectionType::PlanarXY, 1.0f);
    ASSERT(std::fabs(unitQuad.uvs[0] - 0.0f) < 1e-4f && std::fabs(unitQuad.uvs[1] - 0.0f) < 1e-4f,
           "planar XY: (0,0) -> (0,0)");
    ASSERT(std::fabs(unitQuad.uvs[2] - 1.0f) < 1e-4f && std::fabs(unitQuad.uvs[3] - 0.0f) < 1e-4f,
           "planar XY: (1,0) -> (1,0)");
    ASSERT(std::fabs(unitQuad.uvs[4] - 1.0f) < 1e-4f && std::fabs(unitQuad.uvs[5] - 1.0f) < 1e-4f,
           "planar XY: (1,1) -> (1,1)");
    ASSERT(std::fabs(unitQuad.uvs[6] - 0.0f) < 1e-4f && std::fabs(unitQuad.uvs[7] - 1.0f) < 1e-4f,
           "planar XY: (0,1) -> (0,1)");
}

TEST(project_uvs_planar_xz_and_yz_analytic) {
    // Quad on XZ plane
    bromesh::MeshData xzQuad;
    xzQuad.positions = {
        0.0f, 0.0f, 0.0f,
        1.0f, 0.0f, 0.0f,
        1.0f, 0.0f, 1.0f,
        0.0f, 0.0f, 1.0f
    };
    xzQuad.indices = { 0, 1, 2, 0, 2, 3 };
    bromesh::projectUVs(xzQuad, bromesh::ProjectionType::PlanarXZ, 1.0f);
    ASSERT(std::fabs(xzQuad.uvs[0] - 0.0f) < 1e-4f && std::fabs(xzQuad.uvs[1] - 0.0f) < 1e-4f,
           "planar XZ: (0,0,0) -> (0,0)");
    ASSERT(std::fabs(xzQuad.uvs[2] - 1.0f) < 1e-4f && std::fabs(xzQuad.uvs[3] - 0.0f) < 1e-4f,
           "planar XZ: (1,0,0) -> (1,0)");
    ASSERT(std::fabs(xzQuad.uvs[4] - 1.0f) < 1e-4f && std::fabs(xzQuad.uvs[5] - 1.0f) < 1e-4f,
           "planar XZ: (1,0,1) -> (1,1)");
    ASSERT(std::fabs(xzQuad.uvs[6] - 0.0f) < 1e-4f && std::fabs(xzQuad.uvs[7] - 1.0f) < 1e-4f,
           "planar XZ: (0,0,1) -> (0,1)");

    // Quad on YZ plane
    bromesh::MeshData yzQuad;
    yzQuad.positions = {
        0.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f,
        0.0f, 1.0f, 1.0f,
        0.0f, 0.0f, 1.0f
    };
    yzQuad.indices = { 0, 1, 2, 0, 2, 3 };
    bromesh::projectUVs(yzQuad, bromesh::ProjectionType::PlanarYZ, 1.0f);
    ASSERT(std::fabs(yzQuad.uvs[0] - 0.0f) < 1e-4f && std::fabs(yzQuad.uvs[1] - 0.0f) < 1e-4f,
           "planar YZ: (0,0,0) -> (0,0)");
    ASSERT(std::fabs(yzQuad.uvs[2] - 1.0f) < 1e-4f && std::fabs(yzQuad.uvs[3] - 0.0f) < 1e-4f,
           "planar YZ: (0,1,0) -> (1,0)");
    ASSERT(std::fabs(yzQuad.uvs[4] - 1.0f) < 1e-4f && std::fabs(yzQuad.uvs[5] - 1.0f) < 1e-4f,
           "planar YZ: (0,1,1) -> (1,1)");
    ASSERT(std::fabs(yzQuad.uvs[6] - 0.0f) < 1e-4f && std::fabs(yzQuad.uvs[7] - 1.0f) < 1e-4f,
           "planar YZ: (0,0,1) -> (0,1)");
}

TEST(project_uvs_spherical_analytic) {
    const float R = 2.0f;
    const float PI = 3.14159265358979323846f;
    bromesh::MeshData mesh;
    // Points on sphere: (R, 0, 0), (0, 0, R), (-R, 0, 0), (0, 0, -R), (0, R, 0), (0, -R, 0)
    mesh.positions = {
         R,    0.0f,  0.0f, // 0: +X equator
         0.0f, 0.0f,  R,    // 1: +Z equator
        -R,    0.0f,  0.0f, // 2: -X equator
         0.0f, 0.0f, -R,    // 3: -Z equator
         0.0f, R,     0.0f, // 4: +Y North pole
         0.0f,-R,     0.0f  // 5: -Y South pole
    };
    mesh.indices = { 0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0, 4, 0, 5, 1, 1, 5, 2, 2, 5, 3, 3, 5, 0 };

    bromesh::projectUVs(mesh, bromesh::ProjectionType::Spherical, 1.0f);

    // u = atan2(z, x) / (2*PI)
    // v = acos(y / R) / PI
    // Vertex 0: (R, 0, 0) -> u = 0, v = 0.5
    ASSERT(std::fabs(mesh.uvs[0*2+0] - 0.0f) < 1e-4f, "spherical: +X equator u == 0");
    ASSERT(std::fabs(mesh.uvs[0*2+1] - 0.5f) < 1e-4f, "spherical: +X equator v == 0.5");

    // Vertex 1: (0, 0, R) -> u = atan2(R, 0)/(2*PI) = 0.25, v = 0.5
    ASSERT(std::fabs(mesh.uvs[1*2+0] - 0.25f) < 1e-4f, "spherical: +Z equator u == 0.25");
    ASSERT(std::fabs(mesh.uvs[1*2+1] - 0.5f) < 1e-4f, "spherical: +Z equator v == 0.5");

    // Vertex 2: (-R, 0, 0) -> u = atan2(0, -R)/(2*PI) = 0.5, v = 0.5
    ASSERT(std::fabs(mesh.uvs[2*2+0] - 0.5f) < 1e-4f, "spherical: -X equator u == 0.5");
    ASSERT(std::fabs(mesh.uvs[2*2+1] - 0.5f) < 1e-4f, "spherical: -X equator v == 0.5");

    // Vertex 3: (0, 0, -R) -> u = atan2(-R, 0)/(2*PI) = -0.25, v = 0.5
    ASSERT(std::fabs(mesh.uvs[3*2+0] - (-0.25f)) < 1e-4f, "spherical: -Z equator u == -0.25");
    ASSERT(std::fabs(mesh.uvs[3*2+1] - 0.5f) < 1e-4f, "spherical: -Z equator v == 0.5");

    // Vertex 4: (0, R, 0) -> v = acos(1)/PI = 0.0
    ASSERT(std::fabs(mesh.uvs[4*2+1] - 0.0f) < 1e-4f, "spherical: North pole v == 0.0");

    // Vertex 5: (0, -R, 0) -> v = acos(-1)/PI = 1.0
    ASSERT(std::fabs(mesh.uvs[5*2+1] - 1.0f) < 1e-4f, "spherical: South pole v == 1.0");
}

TEST(project_uvs_cylindrical_analytic) {
    const float R = 2.0f;
    bromesh::MeshData mesh;
    mesh.positions = {
         R,    1.5f,  0.0f, // 0: +X, y=1.5
         0.0f, 1.5f,  R,    // 1: +Z, y=1.5
        -R,    0.5f,  0.0f, // 2: -X, y=0.5
         0.0f, 0.5f, -R     // 3: -Z, y=0.5
    };
    mesh.indices = { 0, 1, 2, 0, 2, 3 };

    bromesh::projectUVs(mesh, bromesh::ProjectionType::Cylindrical, 1.0f);
    // u = atan2(z, x)/(2*PI), v = y * scale
    ASSERT(std::fabs(mesh.uvs[0*2+0] - 0.0f) < 1e-4f, "cylindrical: +X u == 0");
    ASSERT(std::fabs(mesh.uvs[0*2+1] - 1.5f) < 1e-4f, "cylindrical: y=1.5 v == 1.5");

    ASSERT(std::fabs(mesh.uvs[1*2+0] - 0.25f) < 1e-4f, "cylindrical: +Z u == 0.25");
    ASSERT(std::fabs(mesh.uvs[1*2+1] - 1.5f) < 1e-4f, "cylindrical: y=1.5 v == 1.5");

    ASSERT(std::fabs(mesh.uvs[2*2+0] - 0.5f) < 1e-4f, "cylindrical: -X u == 0.5");
    ASSERT(std::fabs(mesh.uvs[2*2+1] - 0.5f) < 1e-4f, "cylindrical: y=0.5 v == 0.5");

    ASSERT(std::fabs(mesh.uvs[3*2+0] - (-0.25f)) < 1e-4f, "cylindrical: -Z u == -0.25");
    ASSERT(std::fabs(mesh.uvs[3*2+1] - 0.5f) < 1e-4f, "cylindrical: y=0.5 v == 0.5");
}

TEST(project_uvs_box_analytic) {
    auto b = bromesh::box(1.0f, 1.0f, 1.0f);
    b.uvs.clear();
    bromesh::projectUVs(b, bromesh::ProjectionType::Box, 1.0f);
    ASSERT(b.hasUVs(), "box projection produces UVs");

    // Face 0: +Z face (front) -> dominant Z -> u = x * scale, v = y * scale
    // Box vertices for +Z face are in [-1, 1]
    for (size_t v = 0; v < b.vertexCount(); ++v) {
        float x = b.positions[v*3+0];
        float y = b.positions[v*3+1];
        float z = b.positions[v*3+2];
        float u = b.uvs[v*2+0];
        float w = b.uvs[v*2+1];
        // Coordinates should be in [-1, 1]
        ASSERT(std::fabs(u) <= 1.0001f && std::fabs(w) <= 1.0001f,
               "box projection UV coordinates within [-1, 1]");
    }
}


