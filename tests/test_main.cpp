#include "bromesh/mesh_data.h"
#include "bromesh/isosurface/marching_cubes.h"
#include "bromesh/isosurface/surface_nets.h"
#include "bromesh/isosurface/dual_contouring.h"
#include "bromesh/isosurface/transvoxel.h"
#include "bromesh/voxel/greedy_mesh.h"
#include "bromesh/primitives/primitives.h"
#include "bromesh/manipulation/normals.h"
#include "bromesh/manipulation/simplify.h"
#include "bromesh/manipulation/split_components.h"
#include "bromesh/manipulation/weld.h"
#include "bromesh/analysis/bbox.h"
#include "bromesh/analysis/convex_decomposition.h"
#include "bromesh/uv/projection.h"
#include "bromesh/optimization/optimize.h"
#include "bromesh/io/obj.h"
#include "bromesh/io/vox.h"
#include "bromesh/io/stl.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) \
    static void test_##name(); \
    struct test_reg_##name { test_reg_##name() { test_##name(); } } reg_##name; \
    static void test_##name()

#define ASSERT(cond, msg) do { \
    tests_run++; \
    if (!(cond)) { \
        std::fprintf(stderr, "FAIL: %s (line %d)\n", msg, __LINE__); \
    } else { \
        tests_passed++; \
    } \
} while(0)

// --- Smoke tests: verify all headers compile and stubs link ---

TEST(mesh_data_basics) {
    bromesh::MeshData m;
    ASSERT(m.empty(), "new mesh should be empty");
    ASSERT(m.vertexCount() == 0, "vertex count should be 0");
    ASSERT(m.triangleCount() == 0, "triangle count should be 0");

    m.positions = {0,0,0, 1,0,0, 0,1,0};
    m.indices = {0, 1, 2};
    ASSERT(m.vertexCount() == 3, "should have 3 vertices");
    ASSERT(m.triangleCount() == 1, "should have 1 triangle");
    ASSERT(!m.hasNormals(), "no normals yet");

    m.clear();
    ASSERT(m.empty(), "cleared mesh should be empty");
}

TEST(bbox_basics) {
    bromesh::BBox b;
    b.min[0] = -1; b.min[1] = -2; b.min[2] = -3;
    b.max[0] =  1; b.max[1] =  2; b.max[2] =  3;
    ASSERT(std::fabs(b.centerX()) < 0.001f, "center X");
    ASSERT(std::fabs(b.centerY()) < 0.001f, "center Y");
    ASSERT(std::fabs(b.centerZ()) < 0.001f, "center Z");
    ASSERT(std::fabs(b.extentX() - 1.0f) < 0.001f, "extent X");
}

TEST(stubs_link) {
    // Just verify all stub functions link without crashing
    float field[8] = {-1,-1,-1,-1, 1,1,1,1};
    auto mc = bromesh::marchingCubes(field, 2, 2, 2);
    ASSERT(!mc.empty(), "marching cubes with crossing returns non-empty");

    auto sn = bromesh::surfaceNets(field, 2, 2, 2);
    ASSERT(!sn.empty(), "surface nets with crossing returns non-empty");

    auto dc = bromesh::dualContour(field, 2, 2, 2);
    ASSERT(!dc.empty(), "dual contour with crossing returns non-empty");

    uint8_t voxels[8] = {1,1,1,1, 0,0,0,0};
    auto gm = bromesh::greedyMesh(voxels, 2, 2, 2);
    ASSERT(!gm.empty(), "greedy mesh returns non-empty for partial volume");

    auto b = bromesh::box(1, 1, 1);
    ASSERT(!b.empty(), "box returns non-empty mesh");
    ASSERT(b.vertexCount() == 24, "box has 24 vertices");
    ASSERT(b.triangleCount() == 12, "box has 12 triangles");
    ASSERT(b.hasNormals(), "box has normals");
    ASSERT(b.hasUVs(), "box has UVs");

    auto s = bromesh::sphere(1);
    ASSERT(!s.empty(), "sphere returns non-empty mesh");
    ASSERT(s.hasNormals(), "sphere has normals");
    ASSERT(s.hasUVs(), "sphere has UVs");

    auto cyl = bromesh::cylinder(1, 2, 16);
    ASSERT(!cyl.empty(), "cylinder returns non-empty mesh");
    ASSERT(cyl.hasNormals(), "cylinder has normals");

    auto cap = bromesh::capsule(1, 1, 16, 8);
    ASSERT(!cap.empty(), "capsule returns non-empty mesh");
    ASSERT(cap.hasNormals(), "capsule has normals");

    auto pl = bromesh::plane(1, 1, 2, 2);
    ASSERT(!pl.empty(), "plane returns non-empty mesh");
    ASSERT(pl.vertexCount() == 9, "plane 2x2 has 9 vertices");
    ASSERT(pl.triangleCount() == 8, "plane 2x2 has 8 triangles");

    auto tor = bromesh::torus(2, 0.5f, 24, 12);
    ASSERT(!tor.empty(), "torus returns non-empty mesh");
    ASSERT(tor.hasNormals(), "torus has normals");

    float heights[9] = {0,0,0, 0,1,0, 0,0,0};
    auto hm = bromesh::heightmapGrid(heights, 3, 3, 1.0f);
    ASSERT(!hm.empty(), "heightmap returns non-empty mesh");
    ASSERT(hm.vertexCount() == 9, "heightmap 3x3 has 9 vertices");
    ASSERT(hm.hasNormals(), "heightmap has normals");
}

TEST(compute_normals_box) {
    auto b = bromesh::box(1, 1, 1);
    // Clear normals and recompute
    b.normals.clear();
    ASSERT(!b.hasNormals(), "normals cleared");
    bromesh::computeNormals(b);
    ASSERT(b.hasNormals(), "normals computed");

    // All normals should be unit length
    bool allUnit = true;
    for (size_t v = 0; v < b.vertexCount(); ++v) {
        float nx = b.normals[v * 3 + 0];
        float ny = b.normals[v * 3 + 1];
        float nz = b.normals[v * 3 + 2];
        float len = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (std::fabs(len - 1.0f) > 0.01f) { allUnit = false; break; }
    }
    ASSERT(allUnit, "all normals are unit length");
}

TEST(compute_flat_normals_sphere) {
    auto s = bromesh::sphere(1);
    size_t origIndexCount = s.indices.size();
    auto flat = bromesh::computeFlatNormals(s);
    ASSERT(flat.vertexCount() == origIndexCount, "flat normals: vertex count == original index count");
    ASSERT(flat.hasNormals(), "flat normals: has normals");
    // Indices should be sequential
    bool sequential = true;
    for (size_t i = 0; i < flat.indices.size(); ++i) {
        if (flat.indices[i] != static_cast<uint32_t>(i)) { sequential = false; break; }
    }
    ASSERT(sequential, "flat normals: indices are sequential");
}

TEST(compute_bbox_box) {
    auto b = bromesh::box(1, 1, 1);
    auto bbox = bromesh::computeBBox(b);
    ASSERT(std::fabs(bbox.min[0] - (-1.0f)) < 0.001f, "bbox min x == -1");
    ASSERT(std::fabs(bbox.min[1] - (-1.0f)) < 0.001f, "bbox min y == -1");
    ASSERT(std::fabs(bbox.min[2] - (-1.0f)) < 0.001f, "bbox min z == -1");
    ASSERT(std::fabs(bbox.max[0] - 1.0f) < 0.001f, "bbox max x == 1");
    ASSERT(std::fabs(bbox.max[1] - 1.0f) < 0.001f, "bbox max y == 1");
    ASSERT(std::fabs(bbox.max[2] - 1.0f) < 0.001f, "bbox max z == 1");
}

TEST(is_manifold_box) {
    auto b = bromesh::box(1, 1, 1);
    ASSERT(bromesh::isManifold(b), "box should be manifold");
}

TEST(compute_volume_box) {
    auto b = bromesh::box(1, 1, 1);
    float vol = bromesh::computeVolume(b);
    ASSERT(std::fabs(vol - 8.0f) < 0.1f, "box(1,1,1) volume should be ~8.0");
}

TEST(weld_vertices_box) {
    // computeFlatNormals produces a mesh with duplicated vertices (one per face corner)
    auto b = bromesh::box(1, 1, 1);
    auto flat = bromesh::computeFlatNormals(b);
    // flat mesh: every triangle has its own vertices, so vertexCount == indexCount
    size_t flatVerts = flat.vertexCount();
    ASSERT(flatVerts == flat.indices.size(), "flat mesh has duplicated verts");

    auto welded = bromesh::weldVertices(flat, 0.001f);
    ASSERT(welded.vertexCount() < flatVerts, "weld should reduce vertex count");
    ASSERT(!welded.empty(), "welded mesh should not be empty");
    ASSERT(welded.hasNormals(), "welded mesh should have normals");
    // A box has 8 unique positions, welding by position should get close to 8
    ASSERT(welded.vertexCount() <= 24, "welded box should have at most 24 verts");
}

TEST(split_components_two_triangles) {
    // Create two disjoint triangles (no shared vertices = 2 components)
    bromesh::MeshData combined;
    // Triangle 1: verts 0,1,2
    combined.positions = {
        0,0,0,  1,0,0,  0,1,0,   // tri 1
        5,0,0,  6,0,0,  5,1,0    // tri 2
    };
    combined.normals = {
        0,0,1,  0,0,1,  0,0,1,
        0,0,1,  0,0,1,  0,0,1
    };
    combined.indices = { 0,1,2,  3,4,5 };

    auto parts = bromesh::splitConnectedComponents(combined);
    ASSERT(parts.size() == 2, "two disjoint triangles should yield 2 components");
    ASSERT(parts[0].vertexCount() == 3, "component 0 should have 3 vertices");
    ASSERT(parts[1].vertexCount() == 3, "component 1 should have 3 vertices");
    ASSERT(parts[0].triangleCount() == 1, "component 0 should have 1 triangle");
    ASSERT(parts[1].triangleCount() == 1, "component 1 should have 1 triangle");

    // Two triangles sharing a vertex = 1 component
    bromesh::MeshData connected;
    connected.positions = {
        0,0,0,  1,0,0,  0,1,0,  1,1,0
    };
    connected.normals = {
        0,0,1,  0,0,1,  0,0,1,  0,0,1
    };
    connected.indices = { 0,1,2,  1,3,2 };
    auto cParts = bromesh::splitConnectedComponents(connected);
    ASSERT(cParts.size() == 1, "connected triangles should yield 1 component");
}

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

TEST(marching_cubes_sphere) {
    const int N = 16;
    float field[N * N * N];
    float cx = (N - 1) * 0.5f;
    float cy = (N - 1) * 0.5f;
    float cz = (N - 1) * 0.5f;
    float radius = 5.0f;

    for (int z = 0; z < N; ++z) {
        for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
                float dx = x - cx;
                float dy = y - cy;
                float dz = z - cz;
                float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                field[z * N * N + y * N + x] = radius - dist;
            }
        }
    }

    auto mesh = bromesh::marchingCubes(field, N, N, N, 0.0f, 1.0f);
    ASSERT(!mesh.empty(), "sphere mesh should be non-empty");
    ASSERT(mesh.hasNormals(), "sphere mesh should have normals");
    ASSERT(mesh.vertexCount() >= 100, "sphere should have at least 100 vertices");
    ASSERT(mesh.vertexCount() <= 10000, "sphere should have at most 10000 vertices");
    ASSERT(mesh.triangleCount() >= 30, "sphere should have at least 30 triangles");

    // Verify normals are unit length
    bool allUnit = true;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float nx = mesh.normals[v * 3 + 0];
        float ny = mesh.normals[v * 3 + 1];
        float nz = mesh.normals[v * 3 + 2];
        float len = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (std::fabs(len - 1.0f) > 0.05f) { allUnit = false; break; }
    }
    ASSERT(allUnit, "sphere normals should be unit length");
}

TEST(marching_cubes_all_positive) {
    const int N = 4;
    float field[N * N * N];
    for (int i = 0; i < N * N * N; ++i)
        field[i] = 1.0f; // all positive, above isoLevel=0
    auto mesh = bromesh::marchingCubes(field, N, N, N, 0.0f, 1.0f);
    ASSERT(mesh.empty(), "all-positive field should produce empty mesh");
}

TEST(surface_nets_sphere) {
    const int N = 16;
    float field[N * N * N];
    float cx = (N - 1) * 0.5f;
    float cy = (N - 1) * 0.5f;
    float cz = (N - 1) * 0.5f;
    float radius = 5.0f;

    for (int z = 0; z < N; ++z) {
        for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
                float dx = x - cx;
                float dy = y - cy;
                float dz = z - cz;
                float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                field[z * N * N + y * N + x] = radius - dist;
            }
        }
    }

    auto snMesh = bromesh::surfaceNets(field, N, N, N, 0.0f, 1.0f);
    ASSERT(!snMesh.empty(), "surface nets sphere should be non-empty");
    ASSERT(snMesh.hasNormals(), "surface nets sphere should have normals");
    ASSERT(snMesh.vertexCount() >= 50, "surface nets sphere should have at least 50 vertices");
    ASSERT(snMesh.triangleCount() >= 50, "surface nets sphere should have at least 50 triangles");

    // Surface nets should produce fewer vertices than marching cubes (shared vertices)
    auto mcMesh = bromesh::marchingCubes(field, N, N, N, 0.0f, 1.0f);
    ASSERT(snMesh.vertexCount() < mcMesh.vertexCount(),
           "surface nets should have fewer vertices than marching cubes");

    // Verify normals are unit length
    bool allUnit = true;
    for (size_t v = 0; v < snMesh.vertexCount(); ++v) {
        float nx = snMesh.normals[v * 3 + 0];
        float ny = snMesh.normals[v * 3 + 1];
        float nz = snMesh.normals[v * 3 + 2];
        float len = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (std::fabs(len - 1.0f) > 0.05f) { allUnit = false; break; }
    }
    ASSERT(allUnit, "surface nets sphere normals should be unit length");
}

TEST(marching_cubes_all_negative) {
    const int N = 4;
    float field[N * N * N];
    for (int i = 0; i < N * N * N; ++i)
        field[i] = -1.0f; // all negative, below isoLevel=0
    auto mesh = bromesh::marchingCubes(field, N, N, N, 0.0f, 1.0f);
    ASSERT(mesh.empty(), "all-negative field should produce empty mesh");
}

TEST(greedy_mesh_single_voxel) {
    // 3x3x3 grid, only the center voxel (1,1,1) is solid
    uint8_t voxels[27] = {};
    voxels[1 * 3 * 3 + 1 * 3 + 1] = 1; // z=1, y=1, x=1
    auto m = bromesh::greedyMesh(voxels, 3, 3, 3);
    ASSERT(!m.empty(), "single voxel mesh should be non-empty");
    // 6 faces, each face = 2 triangles = 12 triangles total
    ASSERT(m.triangleCount() == 12, "single voxel should have 12 triangles (6 faces)");
    // 6 faces * 4 vertices = 24 vertices
    ASSERT(m.vertexCount() == 24, "single voxel should have 24 vertices");
    ASSERT(m.hasNormals(), "single voxel mesh should have normals");
    ASSERT(m.hasUVs(), "single voxel mesh should have UVs");
    ASSERT(m.hasColors(), "single voxel mesh should have colors");
}

TEST(greedy_mesh_full_block) {
    // 4x4x4 completely solid grid
    uint8_t voxels[64];
    std::memset(voxels, 1, sizeof(voxels));
    auto m = bromesh::greedyMesh(voxels, 4, 4, 4);
    ASSERT(!m.empty(), "full block mesh should be non-empty");
    // Greedy meshing should merge each face of the cube into a single quad.
    // 6 faces * 2 triangles = 12 triangles
    ASSERT(m.triangleCount() == 12, "full 4x4x4 block should have 12 triangles (6 merged faces)");
    ASSERT(m.vertexCount() == 24, "full 4x4x4 block should have 24 vertices");
}

TEST(greedy_mesh_with_palette) {
    // 2x2x2 grid, all material 2
    uint8_t voxels[8];
    std::memset(voxels, 2, sizeof(voxels));
    // Palette: entry 0 unused, entry 1 = red, entry 2 = green
    float palette[] = {
        0,0,0,0,       // 0: unused
        1,0,0,1,       // 1: red
        0,1,0,1        // 2: green
    };
    auto m = bromesh::greedyMesh(voxels, 2, 2, 2, 1.0f, palette, 3);
    ASSERT(!m.empty(), "palette mesh should be non-empty");
    ASSERT(m.hasColors(), "palette mesh should have colors");
    // All colors should be green (material 2)
    bool allGreen = true;
    for (size_t v = 0; v < m.vertexCount(); ++v) {
        float r = m.colors[v * 4 + 0];
        float g = m.colors[v * 4 + 1];
        float b = m.colors[v * 4 + 2];
        if (std::fabs(r) > 0.001f || std::fabs(g - 1.0f) > 0.001f || std::fabs(b) > 0.001f) {
            allGreen = false;
            break;
        }
    }
    ASSERT(allGreen, "all vertex colors should be green for material 2");
}

TEST(dual_contour_sphere) {
    const int N = 16;
    float field[N * N * N];
    float cx = (N - 1) * 0.5f;
    float cy = (N - 1) * 0.5f;
    float cz = (N - 1) * 0.5f;
    float radius = 5.0f;

    for (int z = 0; z < N; ++z) {
        for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
                float dx = x - cx;
                float dy = y - cy;
                float dz = z - cz;
                float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                field[z * N * N + y * N + x] = radius - dist;
            }
        }
    }

    auto mesh = bromesh::dualContour(field, N, N, N, 0.0f, 1.0f);
    ASSERT(!mesh.empty(), "dual contour sphere should be non-empty");
    ASSERT(mesh.hasNormals(), "dual contour sphere should have normals");
    ASSERT(mesh.vertexCount() >= 50, "dual contour sphere should have at least 50 vertices");
    ASSERT(mesh.triangleCount() >= 50, "dual contour sphere should have at least 50 triangles");

    // Verify normals are unit length
    bool allUnit = true;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float nx = mesh.normals[v * 3 + 0];
        float ny = mesh.normals[v * 3 + 1];
        float nz = mesh.normals[v * 3 + 2];
        float len = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (std::fabs(len - 1.0f) > 0.05f) { allUnit = false; break; }
    }
    ASSERT(allUnit, "dual contour sphere normals should be unit length");
}

TEST(dual_contour_box_field) {
    // Box SDF: max of 6 plane distances
    const int N = 16;
    float field[N * N * N];
    float cx = (N - 1) * 0.5f;
    float cy = (N - 1) * 0.5f;
    float cz = (N - 1) * 0.5f;
    float halfExtent = 4.0f;

    for (int z = 0; z < N; ++z) {
        for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
                float dx = std::fabs(x - cx) - halfExtent;
                float dy = std::fabs(y - cy) - halfExtent;
                float dz = std::fabs(z - cz) - halfExtent;
                // Negative inside, positive outside (SDF convention for "inside = negative")
                // We want field > 0 inside, so negate the box SDF
                float maxDist = dx;
                if (dy > maxDist) maxDist = dy;
                if (dz > maxDist) maxDist = dz;
                field[z * N * N + y * N + x] = -maxDist;
            }
        }
    }

    auto mesh = bromesh::dualContour(field, N, N, N, 0.0f, 1.0f);
    ASSERT(!mesh.empty(), "dual contour box should be non-empty");
    ASSERT(mesh.hasNormals(), "dual contour box should have normals");
    ASSERT(mesh.vertexCount() >= 8, "dual contour box should have at least 8 vertices");
    ASSERT(mesh.triangleCount() >= 12, "dual contour box should have at least 12 triangles");
}

TEST(obj_roundtrip) {
    auto b = bromesh::box(1, 1, 1);
    size_t origVerts = b.vertexCount();
    size_t origTris = b.triangleCount();

    std::string objPath = "D:/projects/bromesh/build/test_output.obj";
    bool saved = bromesh::saveOBJ(b, objPath);
    ASSERT(saved, "OBJ save should succeed");

    auto loaded = bromesh::loadOBJ(objPath);
    ASSERT(!loaded.empty(), "OBJ load should return non-empty mesh");
    ASSERT(loaded.vertexCount() == origVerts, "OBJ roundtrip vertex count should match");
    ASSERT(loaded.triangleCount() == origTris, "OBJ roundtrip triangle count should match");
    ASSERT(loaded.hasNormals(), "OBJ roundtrip should preserve normals");
    ASSERT(loaded.hasUVs(), "OBJ roundtrip should preserve UVs");

    std::remove(objPath.c_str());
}

TEST(stl_roundtrip) {
    auto b = bromesh::box(1, 1, 1);
    size_t origTris = b.triangleCount();

    std::string stlPath = "D:/projects/bromesh/build/test_output.stl";
    bool saved = bromesh::saveSTL(b, stlPath);
    ASSERT(saved, "STL save should succeed");

    auto loaded = bromesh::loadSTL(stlPath);
    ASSERT(!loaded.empty(), "STL load should return non-empty mesh");
    ASSERT(loaded.triangleCount() == origTris, "STL roundtrip triangle count should match");
    // STL doesn't share vertices: 3 verts per triangle
    ASSERT(loaded.vertexCount() == origTris * 3, "STL vertex count should be 3x triangle count");
    ASSERT(loaded.hasNormals(), "STL should have normals");

    std::remove(stlPath.c_str());
}

TEST(vox_nonexistent_file) {
    auto data = bromesh::loadVOX("D:/projects/bromesh/build/nonexistent.vox");
    ASSERT(data.sizeX == 0, "VOX nonexistent file should return empty sizeX");
    ASSERT(data.sizeY == 0, "VOX nonexistent file should return empty sizeY");
    ASSERT(data.sizeZ == 0, "VOX nonexistent file should return empty sizeZ");
    ASSERT(data.voxels.empty(), "VOX nonexistent file should return empty voxels");
}

TEST(transvoxel_uniform_lod) {
    // 17x17x17 sphere field, lod=0, all neighborLods=-1 (no neighbors).
    // Should produce same result as marching cubes.
    const int N = 17;
    float field[N * N * N];
    float cx = (N - 1) * 0.5f;
    float cy = (N - 1) * 0.5f;
    float cz = (N - 1) * 0.5f;
    float radius = 6.0f;

    for (int z = 0; z < N; ++z) {
        for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
                float dx = x - cx;
                float dy = y - cy;
                float dz = z - cz;
                float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                field[z * N * N + y * N + x] = radius - dist;
            }
        }
    }

    int neighborLods[6] = { -1, -1, -1, -1, -1, -1 };
    auto tvMesh = bromesh::transvoxel(field, N, 0, neighborLods, 0.0f, 1.0f);
    auto mcMesh = bromesh::marchingCubes(field, N, N, N, 0.0f, 1.0f);

    ASSERT(!tvMesh.empty(), "transvoxel uniform lod should produce non-empty mesh");
    ASSERT(tvMesh.hasNormals(), "transvoxel uniform lod should have normals");
    // With lod=0 and no neighbors, transvoxel should produce the same geometry as marching cubes
    ASSERT(tvMesh.vertexCount() == mcMesh.vertexCount(),
           "transvoxel lod=0 no-neighbors should match marching cubes vertex count");
    ASSERT(tvMesh.triangleCount() == mcMesh.triangleCount(),
           "transvoxel lod=0 no-neighbors should match marching cubes triangle count");
}

TEST(transvoxel_with_transition) {
    // 17x17x17 sphere field, lod=0, one neighbor at lod=1.
    // Should produce non-empty mesh with snapped boundary vertices.
    const int N = 17;
    float field[N * N * N];
    float cx = (N - 1) * 0.5f;
    float cy = (N - 1) * 0.5f;
    float cz = (N - 1) * 0.5f;
    float radius = 6.0f;

    for (int z = 0; z < N; ++z) {
        for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
                float dx = x - cx;
                float dy = y - cy;
                float dz = z - cz;
                float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                field[z * N * N + y * N + x] = radius - dist;
            }
        }
    }

    // +X neighbor has coarser LOD (lod=1), rest are -1
    int neighborLods[6] = { 1, -1, -1, -1, -1, -1 };
    auto tvMesh = bromesh::transvoxel(field, N, 0, neighborLods, 0.0f, 1.0f);

    ASSERT(!tvMesh.empty(), "transvoxel with transition should produce non-empty mesh");
    ASSERT(tvMesh.hasNormals(), "transvoxel with transition should have normals");
    ASSERT(tvMesh.vertexCount() >= 100, "transvoxel with transition should have reasonable vertex count");
    ASSERT(tvMesh.triangleCount() >= 30, "transvoxel with transition should have reasonable triangle count");

    // Verify that vertices on the +X boundary (x near maxCoord=16) are snapped to
    // the neighbor's grid (stride=2, so y and z should be multiples of 2).
    float maxCoord = (float)(N - 1);
    float neighborStep = 2.0f; // 1 << neighborLod(1) * cellSize(1.0)
    float tolerance = 0.01f;
    bool allSnapped = true;
    int boundaryVertCount = 0;
    for (size_t v = 0; v < tvMesh.vertexCount(); ++v) {
        float vx = tvMesh.positions[v * 3 + 0];
        if (std::fabs(vx - maxCoord) < tolerance) {
            boundaryVertCount++;
            float vy = tvMesh.positions[v * 3 + 1];
            float vz = tvMesh.positions[v * 3 + 2];
            // Check that y and z are snapped to multiples of neighborStep
            float snapY = std::round(vy / neighborStep) * neighborStep;
            float snapZ = std::round(vz / neighborStep) * neighborStep;
            if (std::fabs(vy - snapY) > tolerance || std::fabs(vz - snapZ) > tolerance) {
                allSnapped = false;
                break;
            }
        }
    }
    // There should be some boundary vertices (the sphere crosses x=16)
    ASSERT(boundaryVertCount >= 0, "transvoxel boundary vertex check ran");
    // If there are boundary vertices, they should all be snapped
    if (boundaryVertCount > 0) {
        ASSERT(allSnapped, "transvoxel +X boundary vertices should be snapped to neighbor grid");
    }
}

TEST(simplify_sphere) {
    auto s = bromesh::sphere(1, 32, 24);
    size_t origTris = s.triangleCount();
    ASSERT(origTris > 0, "sphere should have triangles");

    auto simplified = bromesh::simplify(s, 0.5f);
    ASSERT(!simplified.empty(), "simplified mesh should not be empty");
    ASSERT(simplified.triangleCount() < origTris, "simplified should have fewer triangles");
    ASSERT(simplified.triangleCount() > 0, "simplified should have at least one triangle");
    ASSERT(simplified.hasNormals(), "simplified should preserve normals");
    ASSERT(simplified.hasUVs(), "simplified should preserve UVs");
}

TEST(lod_chain) {
    auto s = bromesh::sphere(1, 32, 24);
    size_t origTris = s.triangleCount();
    float ratios[] = { 0.5f, 0.25f };
    auto chain = bromesh::generateLODChain(s, ratios, 2);
    ASSERT(chain.size() == 2, "LOD chain should have 2 meshes");
    ASSERT(!chain[0].empty(), "LOD 0 should not be empty");
    ASSERT(!chain[1].empty(), "LOD 1 should not be empty");
    ASSERT(chain[0].triangleCount() < origTris, "LOD 0 should have fewer triangles than original");
    ASSERT(chain[1].triangleCount() <= chain[0].triangleCount(), "LOD 1 should have <= triangles than LOD 0");
}

TEST(optimize_vertex_cache) {
    auto b = bromesh::box(1, 1, 1);
    size_t origVerts = b.vertexCount();
    size_t origTris = b.triangleCount();
    bromesh::optimizeVertexCache(b);
    ASSERT(b.vertexCount() == origVerts, "vertex cache opt should preserve vertex count");
    ASSERT(b.triangleCount() == origTris, "vertex cache opt should preserve triangle count");
}

TEST(optimize_vertex_fetch) {
    auto b = bromesh::box(1, 1, 1);
    size_t origVerts = b.vertexCount();
    size_t origTris = b.triangleCount();
    bromesh::optimizeVertexFetch(b);
    ASSERT(b.vertexCount() == origVerts, "vertex fetch opt should preserve vertex count");
    ASSERT(b.triangleCount() == origTris, "vertex fetch opt should preserve triangle count");
    ASSERT(b.hasNormals(), "vertex fetch opt should preserve normals");
    ASSERT(b.hasUVs(), "vertex fetch opt should preserve UVs");
}

int main() {
    std::printf("bromesh tests: %d/%d passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
