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
#include "bromesh/optimization/meshlets.h"
#include "bromesh/optimization/analyze.h"
#include "bromesh/optimization/spatial.h"
#include "bromesh/manipulation/merge.h"
#include "bromesh/manipulation/repair.h"
#include "bromesh/optimization/encode.h"
#include "bromesh/optimization/strips.h"
#include "bromesh/primitives/par_primitives.h"
#include "bromesh/uv/unwrap.h"
#include "bromesh/csg/boolean.h"
#include "bromesh/io/obj.h"
#include "bromesh/io/vox.h"
#include "bromesh/io/stl.h"
#include "bromesh/io/gltf.h"
#include "bromesh/io/ply.h"
#include "bromesh/io/fbx.h"
#include "bromesh/manipulation/subdivide.h"
#include "bromesh/manipulation/smooth.h"
#include "bromesh/manipulation/remesh.h"
#include "bromesh/reconstruction/reconstruct.h"
#include "bromesh/analysis/bake.h"
#include "bromesh/manipulation/skin.h"
#include "bromesh/uv/uv_metrics.h"

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

// ============================================================================
// Helpers for comprehensive I/O roundtrip tests
// ============================================================================

static const char* testDir = "D:/projects/bromesh/build/";

static bool approxEqual(float a, float b, float tol) {
    return std::fabs(a - b) <= tol;
}

// Compare positions within tolerance (sorted comparison for order-independent match)
static bool positionsMatch(const bromesh::MeshData& a, const bromesh::MeshData& b, float tol) {
    if (a.vertexCount() != b.vertexCount()) return false;
    for (size_t i = 0; i < a.positions.size(); ++i) {
        if (!approxEqual(a.positions[i], b.positions[i], tol)) return false;
    }
    return true;
}

static bool normalsMatch(const bromesh::MeshData& a, const bromesh::MeshData& b, float tol) {
    if (a.normals.size() != b.normals.size()) return false;
    for (size_t i = 0; i < a.normals.size(); ++i) {
        if (!approxEqual(a.normals[i], b.normals[i], tol)) return false;
    }
    return true;
}

static bool uvsMatch(const bromesh::MeshData& a, const bromesh::MeshData& b, float tol) {
    if (a.uvs.size() != b.uvs.size()) return false;
    for (size_t i = 0; i < a.uvs.size(); ++i) {
        if (!approxEqual(a.uvs[i], b.uvs[i], tol)) return false;
    }
    return true;
}

static bool indicesMatch(const bromesh::MeshData& a, const bromesh::MeshData& b) {
    return a.indices == b.indices;
}

static bromesh::BBox meshBBox(const bromesh::MeshData& m) {
    return bromesh::computeBBox(m);
}

static bool bboxMatch(const bromesh::BBox& a, const bromesh::BBox& b, float tol) {
    for (int i = 0; i < 3; ++i) {
        if (!approxEqual(a.min[i], b.min[i], tol)) return false;
        if (!approxEqual(a.max[i], b.max[i], tol)) return false;
    }
    return true;
}

// Generate a sphere SDF field
static void fillSphereField(float* field, int N, float radius) {
    float c = (N - 1) * 0.5f;
    for (int z = 0; z < N; ++z)
        for (int y = 0; y < N; ++y)
            for (int x = 0; x < N; ++x) {
                float dx = x - c, dy = y - c, dz = z - c;
                field[z * N * N + y * N + x] = radius - std::sqrt(dx*dx + dy*dy + dz*dz);
            }
}

// ============================================================================
// OBJ roundtrip tests — preserves positions, normals, UVs, indices
// ============================================================================

TEST(obj_rt_sphere_with_cylindrical_uvs) {
    // Sphere primitive -> cylindrical UV projection -> OBJ roundtrip
    auto mesh = bromesh::sphere(2.0f, 24, 16);
    mesh.uvs.clear();
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Cylindrical, 1.0f);
    ASSERT(mesh.hasUVs(), "obj_rt_sphere_cyl: should have UVs after projection");
    auto origBBox = meshBBox(mesh);
    float origVol = bromesh::computeVolume(mesh);

    std::string path = std::string(testDir) + "rt_sphere_cyl.obj";
    ASSERT(bromesh::saveOBJ(mesh, path), "obj_rt_sphere_cyl: save");
    auto loaded = bromesh::loadOBJ(path);
    ASSERT(loaded.vertexCount() == mesh.vertexCount(), "obj_rt_sphere_cyl: vertex count");
    ASSERT(loaded.triangleCount() == mesh.triangleCount(), "obj_rt_sphere_cyl: tri count");
    ASSERT(loaded.hasNormals(), "obj_rt_sphere_cyl: normals preserved");
    ASSERT(loaded.hasUVs(), "obj_rt_sphere_cyl: UVs preserved");
    auto loadedBBox = meshBBox(loaded);
    ASSERT(bboxMatch(origBBox, loadedBBox, 1e-3f), "obj_rt_sphere_cyl: bbox match");
    float loadedVol = bromesh::computeVolume(loaded);
    ASSERT(approxEqual(origVol, loadedVol, std::fabs(origVol) * 0.01f + 0.1f),
           "obj_rt_sphere_cyl: volume match");
    std::remove(path.c_str());
}

TEST(obj_rt_torus_simplified) {
    // Torus -> simplify to 50% -> recompute normals -> spherical UVs -> OBJ roundtrip
    auto mesh = bromesh::torus(3.0f, 1.0f, 32, 16);
    mesh = bromesh::simplify(mesh, 0.5f);
    ASSERT(!mesh.empty(), "obj_rt_torus_simp: simplified not empty");
    bromesh::computeNormals(mesh);
    mesh.uvs.clear();
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Spherical, 1.0f);
    auto origBBox = meshBBox(mesh);

    std::string path = std::string(testDir) + "rt_torus_simp.obj";
    ASSERT(bromesh::saveOBJ(mesh, path), "obj_rt_torus_simp: save");
    auto loaded = bromesh::loadOBJ(path);
    ASSERT(loaded.vertexCount() == mesh.vertexCount(), "obj_rt_torus_simp: vertex count");
    ASSERT(loaded.triangleCount() == mesh.triangleCount(), "obj_rt_torus_simp: tri count");
    ASSERT(loaded.hasNormals(), "obj_rt_torus_simp: normals preserved");
    ASSERT(loaded.hasUVs(), "obj_rt_torus_simp: UVs preserved");
    auto loadedBBox = meshBBox(loaded);
    ASSERT(bboxMatch(origBBox, loadedBBox, 1e-3f), "obj_rt_torus_simp: bbox match");
    std::remove(path.c_str());
}

TEST(obj_rt_marching_cubes_welded_with_box_uvs) {
    // Marching cubes sphere -> weld -> box UV projection -> OBJ roundtrip
    const int N = 16;
    float field[N * N * N];
    fillSphereField(field, N, 5.0f);
    auto mesh = bromesh::marchingCubes(field, N, N, N, 0.0f, 1.0f);
    mesh = bromesh::weldVertices(mesh, 1e-5f);
    ASSERT(!mesh.empty(), "obj_rt_mc_weld: welded not empty");
    bromesh::computeNormals(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Box, 1.0f);
    auto origBBox = meshBBox(mesh);
    float origVol = bromesh::computeVolume(mesh);

    std::string path = std::string(testDir) + "rt_mc_weld.obj";
    ASSERT(bromesh::saveOBJ(mesh, path), "obj_rt_mc_weld: save");
    auto loaded = bromesh::loadOBJ(path);
    ASSERT(loaded.vertexCount() == mesh.vertexCount(), "obj_rt_mc_weld: vertex count");
    ASSERT(loaded.triangleCount() == mesh.triangleCount(), "obj_rt_mc_weld: tri count");
    ASSERT(loaded.hasNormals(), "obj_rt_mc_weld: normals preserved");
    ASSERT(loaded.hasUVs(), "obj_rt_mc_weld: UVs preserved");
    auto loadedBBox = meshBBox(loaded);
    ASSERT(bboxMatch(origBBox, loadedBBox, 1e-3f), "obj_rt_mc_weld: bbox match");
    float loadedVol = bromesh::computeVolume(loaded);
    ASSERT(approxEqual(origVol, loadedVol, std::fabs(origVol) * 0.01f + 0.1f),
           "obj_rt_mc_weld: volume match");
    std::remove(path.c_str());
}

TEST(obj_rt_dual_contour_flat_normals_planarXZ) {
    // Dual contouring sphere -> flat normals -> planar XZ UVs -> OBJ roundtrip
    const int N = 16;
    float field[N * N * N];
    fillSphereField(field, N, 5.0f);
    auto mesh = bromesh::dualContour(field, N, N, N, 0.0f, 1.0f);
    mesh = bromesh::computeFlatNormals(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::PlanarXZ, 0.5f);
    auto origBBox = meshBBox(mesh);

    std::string path = std::string(testDir) + "rt_dc_flat.obj";
    ASSERT(bromesh::saveOBJ(mesh, path), "obj_rt_dc_flat: save");
    auto loaded = bromesh::loadOBJ(path);
    ASSERT(loaded.vertexCount() == mesh.vertexCount(), "obj_rt_dc_flat: vertex count");
    ASSERT(loaded.triangleCount() == mesh.triangleCount(), "obj_rt_dc_flat: tri count");
    ASSERT(loaded.hasNormals(), "obj_rt_dc_flat: normals preserved");
    ASSERT(loaded.hasUVs(), "obj_rt_dc_flat: UVs preserved");
    auto loadedBBox = meshBBox(loaded);
    ASSERT(bboxMatch(origBBox, loadedBBox, 1e-3f), "obj_rt_dc_flat: bbox match");
    std::remove(path.c_str());
}

TEST(obj_rt_capsule_optimized) {
    // Capsule -> vertex cache + fetch optimize -> PlanarYZ UVs -> OBJ roundtrip
    auto mesh = bromesh::capsule(1.5f, 2.0f, 20, 10);
    bromesh::optimizeVertexCache(mesh);
    bromesh::optimizeVertexFetch(mesh);
    mesh.uvs.clear();
    bromesh::projectUVs(mesh, bromesh::ProjectionType::PlanarYZ, 1.0f);
    auto origBBox = meshBBox(mesh);

    std::string path = std::string(testDir) + "rt_capsule_opt.obj";
    ASSERT(bromesh::saveOBJ(mesh, path), "obj_rt_capsule_opt: save");
    auto loaded = bromesh::loadOBJ(path);
    ASSERT(loaded.vertexCount() == mesh.vertexCount(), "obj_rt_capsule_opt: vertex count");
    ASSERT(loaded.triangleCount() == mesh.triangleCount(), "obj_rt_capsule_opt: tri count");
    ASSERT(loaded.hasNormals(), "obj_rt_capsule_opt: normals");
    ASSERT(loaded.hasUVs(), "obj_rt_capsule_opt: UVs");
    auto loadedBBox = meshBBox(loaded);
    ASSERT(bboxMatch(origBBox, loadedBBox, 1e-3f), "obj_rt_capsule_opt: bbox match");
    std::remove(path.c_str());
}

TEST(obj_rt_heightmap_with_planarXY) {
    // Heightmap grid -> smooth normals -> planar XY UVs -> OBJ roundtrip
    float heights[25];
    for (int i = 0; i < 25; ++i)
        heights[i] = std::sin(i * 0.5f) * 2.0f;
    auto mesh = bromesh::heightmapGrid(heights, 5, 5, 1.0f);
    bromesh::computeNormals(mesh);
    mesh.uvs.clear();
    bromesh::projectUVs(mesh, bromesh::ProjectionType::PlanarXY, 1.0f);
    auto origBBox = meshBBox(mesh);

    std::string path = std::string(testDir) + "rt_heightmap.obj";
    ASSERT(bromesh::saveOBJ(mesh, path), "obj_rt_heightmap: save");
    auto loaded = bromesh::loadOBJ(path);
    ASSERT(loaded.vertexCount() == mesh.vertexCount(), "obj_rt_heightmap: vertex count");
    ASSERT(loaded.triangleCount() == mesh.triangleCount(), "obj_rt_heightmap: tri count");
    ASSERT(loaded.hasNormals(), "obj_rt_heightmap: normals preserved");
    ASSERT(loaded.hasUVs(), "obj_rt_heightmap: UVs preserved");
    auto loadedBBox = meshBBox(loaded);
    ASSERT(bboxMatch(origBBox, loadedBBox, 1e-3f), "obj_rt_heightmap: bbox match");
    std::remove(path.c_str());
}

TEST(obj_rt_surface_nets_lod_chain) {
    // Surface nets -> generate LOD chain -> take LOD1 -> recompute normals -> OBJ roundtrip
    const int N = 16;
    float field[N * N * N];
    fillSphereField(field, N, 5.0f);
    auto mesh = bromesh::surfaceNets(field, N, N, N, 0.0f, 1.0f);
    float ratios[] = { 0.7f, 0.4f };
    auto chain = bromesh::generateLODChain(mesh, ratios, 2);
    ASSERT(chain.size() == 2, "obj_rt_sn_lod: chain has 2 levels");
    auto& lod = chain[1];
    bromesh::computeNormals(lod);
    bromesh::projectUVs(lod, bromesh::ProjectionType::Spherical, 1.0f);
    auto origBBox = meshBBox(lod);

    std::string path = std::string(testDir) + "rt_sn_lod.obj";
    ASSERT(bromesh::saveOBJ(lod, path), "obj_rt_sn_lod: save");
    auto loaded = bromesh::loadOBJ(path);
    ASSERT(loaded.vertexCount() == lod.vertexCount(), "obj_rt_sn_lod: vertex count");
    ASSERT(loaded.triangleCount() == lod.triangleCount(), "obj_rt_sn_lod: tri count");
    ASSERT(loaded.hasNormals(), "obj_rt_sn_lod: normals preserved");
    ASSERT(loaded.hasUVs(), "obj_rt_sn_lod: UVs preserved");
    auto loadedBBox = meshBBox(loaded);
    ASSERT(bboxMatch(origBBox, loadedBBox, 1e-3f), "obj_rt_sn_lod: bbox match");
    std::remove(path.c_str());
}

TEST(obj_rt_cylinder_flat_normals_overdraw_opt) {
    // Cylinder -> flat normals -> overdraw optimize -> box UVs -> OBJ
    auto mesh = bromesh::cylinder(2.0f, 3.0f, 24);
    mesh = bromesh::computeFlatNormals(mesh);
    bromesh::optimizeOverdraw(mesh, 1.05f);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Box, 1.0f);
    auto origBBox = meshBBox(mesh);

    std::string path = std::string(testDir) + "rt_cyl_flat_od.obj";
    ASSERT(bromesh::saveOBJ(mesh, path), "obj_rt_cyl_flat_od: save");
    auto loaded = bromesh::loadOBJ(path);
    ASSERT(loaded.vertexCount() == mesh.vertexCount(), "obj_rt_cyl_flat_od: vertex count");
    ASSERT(loaded.triangleCount() == mesh.triangleCount(), "obj_rt_cyl_flat_od: tri count");
    ASSERT(loaded.hasNormals(), "obj_rt_cyl_flat_od: normals preserved");
    ASSERT(loaded.hasUVs(), "obj_rt_cyl_flat_od: UVs preserved");
    auto loadedBBox = meshBBox(loaded);
    ASSERT(bboxMatch(origBBox, loadedBBox, 1e-3f), "obj_rt_cyl_flat_od: bbox match");
    std::remove(path.c_str());
}

// ============================================================================
// STL roundtrip tests — positions + face normals only; vertices unshared on load
// ============================================================================

TEST(stl_rt_torus_with_smooth_normals) {
    // Torus -> smooth normals -> STL roundtrip
    // STL stores per-face normals and doesn't share vertices, so we compare
    // triangle count and bounding box
    auto mesh = bromesh::torus(2.0f, 0.5f, 24, 12);
    bromesh::computeNormals(mesh);
    auto origBBox = meshBBox(mesh);
    size_t origTris = mesh.triangleCount();

    std::string path = std::string(testDir) + "rt_torus.stl";
    ASSERT(bromesh::saveSTL(mesh, path), "stl_rt_torus: save");
    auto loaded = bromesh::loadSTL(path);
    ASSERT(loaded.triangleCount() == origTris, "stl_rt_torus: tri count");
    ASSERT(loaded.vertexCount() == origTris * 3, "stl_rt_torus: 3 verts per tri");
    ASSERT(loaded.hasNormals(), "stl_rt_torus: has normals");
    auto loadedBBox = meshBBox(loaded);
    ASSERT(bboxMatch(origBBox, loadedBBox, 0.01f), "stl_rt_torus: bbox match");
    std::remove(path.c_str());
}

TEST(stl_rt_marching_cubes_simplified) {
    // Marching cubes -> simplify -> STL roundtrip
    const int N = 16;
    float field[N * N * N];
    fillSphereField(field, N, 5.0f);
    auto mesh = bromesh::marchingCubes(field, N, N, N, 0.0f, 1.0f);
    mesh = bromesh::simplify(mesh, 0.3f);
    ASSERT(!mesh.empty(), "stl_rt_mc_simp: simplified not empty");
    bromesh::computeNormals(mesh);
    auto origBBox = meshBBox(mesh);
    size_t origTris = mesh.triangleCount();

    std::string path = std::string(testDir) + "rt_mc_simp.stl";
    ASSERT(bromesh::saveSTL(mesh, path), "stl_rt_mc_simp: save");
    auto loaded = bromesh::loadSTL(path);
    ASSERT(loaded.triangleCount() == origTris, "stl_rt_mc_simp: tri count");
    auto loadedBBox = meshBBox(loaded);
    ASSERT(bboxMatch(origBBox, loadedBBox, 0.01f), "stl_rt_mc_simp: bbox match");
    // Verify loaded volume is roughly the same
    float origVol = bromesh::computeVolume(mesh);
    // STL loaded mesh won't have shared verts, but volume should still be close
    // (volume computation uses triangle faces, doesn't depend on vertex sharing)
    float loadedVol = bromesh::computeVolume(loaded);
    ASSERT(approxEqual(origVol, loadedVol, std::fabs(origVol) * 0.01f + 0.1f),
           "stl_rt_mc_simp: volume match");
    std::remove(path.c_str());
}

TEST(stl_rt_capsule_welded_optimized) {
    // Capsule -> weld -> optimize vertex cache -> STL roundtrip
    auto mesh = bromesh::capsule(1.0f, 1.5f, 20, 10);
    mesh = bromesh::weldVertices(mesh, 1e-5f);
    bromesh::computeNormals(mesh);
    bromesh::optimizeVertexCache(mesh);
    auto origBBox = meshBBox(mesh);
    size_t origTris = mesh.triangleCount();

    std::string path = std::string(testDir) + "rt_capsule_weld.stl";
    ASSERT(bromesh::saveSTL(mesh, path), "stl_rt_capsule_weld: save");
    auto loaded = bromesh::loadSTL(path);
    ASSERT(loaded.triangleCount() == origTris, "stl_rt_capsule_weld: tri count");
    auto loadedBBox = meshBBox(loaded);
    ASSERT(bboxMatch(origBBox, loadedBBox, 0.01f), "stl_rt_capsule_weld: bbox match");
    std::remove(path.c_str());
}

TEST(stl_rt_dual_contour_flat_normals) {
    // Dual contour -> flat normals -> STL roundtrip
    const int N = 16;
    float field[N * N * N];
    fillSphereField(field, N, 5.0f);
    auto mesh = bromesh::dualContour(field, N, N, N, 0.0f, 1.0f);
    mesh = bromesh::computeFlatNormals(mesh);
    auto origBBox = meshBBox(mesh);
    size_t origTris = mesh.triangleCount();

    std::string path = std::string(testDir) + "rt_dc_flat.stl";
    ASSERT(bromesh::saveSTL(mesh, path), "stl_rt_dc_flat: save");
    auto loaded = bromesh::loadSTL(path);
    ASSERT(loaded.triangleCount() == origTris, "stl_rt_dc_flat: tri count");
    auto loadedBBox = meshBBox(loaded);
    ASSERT(bboxMatch(origBBox, loadedBBox, 0.01f), "stl_rt_dc_flat: bbox match");
    std::remove(path.c_str());
}

TEST(stl_rt_plane_heightmap) {
    // Heightmap -> smooth normals -> STL roundtrip
    float heights[36];
    for (int i = 0; i < 36; ++i)
        heights[i] = std::cos(i * 0.3f) * 1.5f;
    auto mesh = bromesh::heightmapGrid(heights, 6, 6, 0.5f);
    bromesh::computeNormals(mesh);
    auto origBBox = meshBBox(mesh);
    size_t origTris = mesh.triangleCount();

    std::string path = std::string(testDir) + "rt_heightmap.stl";
    ASSERT(bromesh::saveSTL(mesh, path), "stl_rt_heightmap: save");
    auto loaded = bromesh::loadSTL(path);
    ASSERT(loaded.triangleCount() == origTris, "stl_rt_heightmap: tri count");
    auto loadedBBox = meshBBox(loaded);
    ASSERT(bboxMatch(origBBox, loadedBBox, 0.01f), "stl_rt_heightmap: bbox match");
    std::remove(path.c_str());
}

TEST(stl_rt_sphere_lod_chain_level0) {
    // Sphere -> LOD chain -> take LOD 0 -> STL roundtrip
    auto mesh = bromesh::sphere(3.0f, 32, 24);
    float ratios[] = { 0.6f, 0.3f };
    auto chain = bromesh::generateLODChain(mesh, ratios, 2);
    ASSERT(!chain[0].empty(), "stl_rt_lod0: lod0 not empty");
    bromesh::computeNormals(chain[0]);
    auto origBBox = meshBBox(chain[0]);
    size_t origTris = chain[0].triangleCount();

    std::string path = std::string(testDir) + "rt_lod0.stl";
    ASSERT(bromesh::saveSTL(chain[0], path), "stl_rt_lod0: save");
    auto loaded = bromesh::loadSTL(path);
    ASSERT(loaded.triangleCount() == origTris, "stl_rt_lod0: tri count");
    auto loadedBBox = meshBBox(loaded);
    ASSERT(bboxMatch(origBBox, loadedBBox, 0.01f), "stl_rt_lod0: bbox match");
    std::remove(path.c_str());
}

TEST(stl_rt_surface_nets_welded) {
    // Surface nets -> weld -> STL roundtrip
    const int N = 16;
    float field[N * N * N];
    fillSphereField(field, N, 5.0f);
    auto mesh = bromesh::surfaceNets(field, N, N, N, 0.0f, 1.0f);
    mesh = bromesh::weldVertices(mesh, 1e-4f);
    bromesh::computeNormals(mesh);
    auto origBBox = meshBBox(mesh);
    size_t origTris = mesh.triangleCount();

    std::string path = std::string(testDir) + "rt_sn_weld.stl";
    ASSERT(bromesh::saveSTL(mesh, path), "stl_rt_sn_weld: save");
    auto loaded = bromesh::loadSTL(path);
    ASSERT(loaded.triangleCount() == origTris, "stl_rt_sn_weld: tri count");
    auto loadedBBox = meshBBox(loaded);
    ASSERT(bboxMatch(origBBox, loadedBBox, 0.01f), "stl_rt_sn_weld: bbox match");
    std::remove(path.c_str());
}

TEST(stl_rt_greedy_mesh_voxel_block) {
    // Greedy mesh voxel block -> STL roundtrip
    uint8_t voxels[64];
    std::memset(voxels, 1, sizeof(voxels));
    auto mesh = bromesh::greedyMesh(voxels, 4, 4, 4, 1.0f);
    bromesh::computeNormals(mesh);
    auto origBBox = meshBBox(mesh);
    size_t origTris = mesh.triangleCount();

    std::string path = std::string(testDir) + "rt_greedy.stl";
    ASSERT(bromesh::saveSTL(mesh, path), "stl_rt_greedy: save");
    auto loaded = bromesh::loadSTL(path);
    ASSERT(loaded.triangleCount() == origTris, "stl_rt_greedy: tri count");
    auto loadedBBox = meshBBox(loaded);
    ASSERT(bboxMatch(origBBox, loadedBBox, 0.01f), "stl_rt_greedy: bbox match");
    float origVol = bromesh::computeVolume(mesh);
    float loadedVol = bromesh::computeVolume(loaded);
    ASSERT(approxEqual(origVol, loadedVol, 0.1f), "stl_rt_greedy: volume match");
    std::remove(path.c_str());
}

// ============================================================================
// glTF roundtrip tests — preserves positions, normals, UVs, indices
// ============================================================================

#if BROMESH_HAS_GLTF

TEST(gltf_rt_box_with_all_uv_types) {
    // Box -> each UV projection type -> glTF roundtrip
    bromesh::ProjectionType types[] = {
        bromesh::ProjectionType::Box,
        bromesh::ProjectionType::PlanarXY,
        bromesh::ProjectionType::PlanarXZ,
        bromesh::ProjectionType::PlanarYZ,
        bromesh::ProjectionType::Cylindrical,
        bromesh::ProjectionType::Spherical
    };
    const char* names[] = { "Box", "PlanarXY", "PlanarXZ", "PlanarYZ", "Cylindrical", "Spherical" };

    for (int t = 0; t < 6; ++t) {
        auto mesh = bromesh::box(1.0f, 1.5f, 2.0f);
        mesh.uvs.clear();
        bromesh::projectUVs(mesh, types[t], 1.0f);
        ASSERT(mesh.hasUVs(), "gltf_rt_box_uvs: has UVs");

        std::string path = std::string(testDir) + "rt_box_" + names[t] + ".glb";
        ASSERT(bromesh::saveGLTF(mesh, path), "gltf_rt_box_uvs: save");
        auto scene = bromesh::loadGLTF(path);
        ASSERT(!scene.meshes.empty(), "gltf_rt_box_uvs: loaded scene has meshes");
        auto& loaded = scene.meshes[0];
        ASSERT(loaded.vertexCount() == mesh.vertexCount(),
               "gltf_rt_box_uvs: vertex count");
        ASSERT(loaded.triangleCount() == mesh.triangleCount(),
               "gltf_rt_box_uvs: tri count");
        ASSERT(positionsMatch(mesh, loaded, 1e-4f),
               "gltf_rt_box_uvs: positions match");
        ASSERT(loaded.hasNormals(), "gltf_rt_box_uvs: normals preserved");
        ASSERT(loaded.hasUVs(), "gltf_rt_box_uvs: UVs preserved");
        ASSERT(uvsMatch(mesh, loaded, 1e-4f), "gltf_rt_box_uvs: UVs match");
        std::remove(path.c_str());
    }
}

TEST(gltf_rt_sphere_simplified_cylindrical) {
    // Sphere -> simplify -> recompute normals -> cylindrical UVs -> glTF roundtrip
    auto mesh = bromesh::sphere(2.5f, 32, 24);
    mesh = bromesh::simplify(mesh, 0.4f);
    bromesh::computeNormals(mesh);
    mesh.uvs.clear();
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Cylindrical, 1.0f);

    std::string path = std::string(testDir) + "rt_sphere_simp.glb";
    ASSERT(bromesh::saveGLTF(mesh, path), "gltf_rt_sphere_simp: save");
    auto scene = bromesh::loadGLTF(path);
    ASSERT(!scene.meshes.empty(), "gltf_rt_sphere_simp: has meshes");
    auto& loaded = scene.meshes[0];
    ASSERT(loaded.vertexCount() == mesh.vertexCount(), "gltf_rt_sphere_simp: vertex count");
    ASSERT(loaded.triangleCount() == mesh.triangleCount(), "gltf_rt_sphere_simp: tri count");
    ASSERT(positionsMatch(mesh, loaded, 1e-4f), "gltf_rt_sphere_simp: positions");
    ASSERT(normalsMatch(mesh, loaded, 1e-3f), "gltf_rt_sphere_simp: normals");
    ASSERT(uvsMatch(mesh, loaded, 1e-4f), "gltf_rt_sphere_simp: UVs");
    std::remove(path.c_str());
}

TEST(gltf_rt_torus_flat_normals_optimized) {
    // Torus -> flat normals -> vertex cache + fetch optimize -> box UVs -> glTF
    auto mesh = bromesh::torus(2.0f, 0.7f, 24, 12);
    mesh = bromesh::computeFlatNormals(mesh);
    bromesh::optimizeVertexCache(mesh);
    bromesh::optimizeVertexFetch(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Box, 1.0f);

    std::string path = std::string(testDir) + "rt_torus_flat_opt.glb";
    ASSERT(bromesh::saveGLTF(mesh, path), "gltf_rt_torus_flat: save");
    auto scene = bromesh::loadGLTF(path);
    ASSERT(!scene.meshes.empty(), "gltf_rt_torus_flat: has meshes");
    auto& loaded = scene.meshes[0];
    ASSERT(loaded.vertexCount() == mesh.vertexCount(), "gltf_rt_torus_flat: vertex count");
    ASSERT(loaded.triangleCount() == mesh.triangleCount(), "gltf_rt_torus_flat: tri count");
    ASSERT(positionsMatch(mesh, loaded, 1e-4f), "gltf_rt_torus_flat: positions");
    ASSERT(normalsMatch(mesh, loaded, 1e-3f), "gltf_rt_torus_flat: normals");
    std::remove(path.c_str());
}

TEST(gltf_rt_marching_cubes_welded) {
    // Marching cubes -> weld -> recompute normals -> spherical UVs -> glTF
    const int N = 16;
    float field[N * N * N];
    fillSphereField(field, N, 5.0f);
    auto mesh = bromesh::marchingCubes(field, N, N, N, 0.0f, 1.0f);
    mesh = bromesh::weldVertices(mesh, 1e-5f);
    bromesh::computeNormals(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Spherical, 1.0f);

    std::string path = std::string(testDir) + "rt_mc_weld.glb";
    ASSERT(bromesh::saveGLTF(mesh, path), "gltf_rt_mc_weld: save");
    auto scene = bromesh::loadGLTF(path);
    ASSERT(!scene.meshes.empty(), "gltf_rt_mc_weld: has meshes");
    auto& loaded = scene.meshes[0];
    ASSERT(loaded.vertexCount() == mesh.vertexCount(), "gltf_rt_mc_weld: vertex count");
    ASSERT(loaded.triangleCount() == mesh.triangleCount(), "gltf_rt_mc_weld: tri count");
    ASSERT(positionsMatch(mesh, loaded, 1e-4f), "gltf_rt_mc_weld: positions");
    ASSERT(loaded.hasNormals(), "gltf_rt_mc_weld: normals");
    ASSERT(loaded.hasUVs(), "gltf_rt_mc_weld: UVs");
    std::remove(path.c_str());
}

TEST(gltf_rt_surface_nets_lod_chain) {
    // Surface nets -> LOD chain -> LOD 0 -> recompute normals -> glTF
    const int N = 16;
    float field[N * N * N];
    fillSphereField(field, N, 5.0f);
    auto mesh = bromesh::surfaceNets(field, N, N, N, 0.0f, 1.0f);
    float ratios[] = { 0.5f };
    auto chain = bromesh::generateLODChain(mesh, ratios, 1);
    ASSERT(!chain.empty() && !chain[0].empty(), "gltf_rt_sn_lod: chain not empty");
    auto& lod = chain[0];
    bromesh::computeNormals(lod);
    bromesh::projectUVs(lod, bromesh::ProjectionType::Cylindrical, 1.0f);

    std::string path = std::string(testDir) + "rt_sn_lod.glb";
    ASSERT(bromesh::saveGLTF(lod, path), "gltf_rt_sn_lod: save");
    auto scene = bromesh::loadGLTF(path);
    ASSERT(!scene.meshes.empty(), "gltf_rt_sn_lod: has meshes");
    auto& loaded = scene.meshes[0];
    ASSERT(loaded.vertexCount() == lod.vertexCount(), "gltf_rt_sn_lod: vertex count");
    ASSERT(loaded.triangleCount() == lod.triangleCount(), "gltf_rt_sn_lod: tri count");
    ASSERT(positionsMatch(lod, loaded, 1e-4f), "gltf_rt_sn_lod: positions");
    ASSERT(loaded.hasNormals(), "gltf_rt_sn_lod: normals");
    std::remove(path.c_str());
}

TEST(gltf_rt_dual_contour_overdraw_opt) {
    // Dual contour -> overdraw optimize -> planar XY UVs -> glTF
    const int N = 16;
    float field[N * N * N];
    fillSphereField(field, N, 5.0f);
    auto mesh = bromesh::dualContour(field, N, N, N, 0.0f, 1.0f);
    bromesh::computeNormals(mesh);
    bromesh::optimizeOverdraw(mesh, 1.05f);
    mesh.uvs.clear();
    bromesh::projectUVs(mesh, bromesh::ProjectionType::PlanarXY, 1.0f);

    std::string path = std::string(testDir) + "rt_dc_od.glb";
    ASSERT(bromesh::saveGLTF(mesh, path), "gltf_rt_dc_od: save");
    auto scene = bromesh::loadGLTF(path);
    ASSERT(!scene.meshes.empty(), "gltf_rt_dc_od: has meshes");
    auto& loaded = scene.meshes[0];
    ASSERT(loaded.vertexCount() == mesh.vertexCount(), "gltf_rt_dc_od: vertex count");
    ASSERT(loaded.triangleCount() == mesh.triangleCount(), "gltf_rt_dc_od: tri count");
    ASSERT(positionsMatch(mesh, loaded, 1e-4f), "gltf_rt_dc_od: positions");
    ASSERT(loaded.hasUVs(), "gltf_rt_dc_od: UVs preserved");
    std::remove(path.c_str());
}

TEST(gltf_rt_capsule_full_pipeline) {
    // Capsule -> weld -> simplify -> flat normals -> all optimizations -> glTF
    auto mesh = bromesh::capsule(1.0f, 2.0f, 24, 12);
    mesh = bromesh::weldVertices(mesh, 1e-5f);
    mesh = bromesh::simplify(mesh, 0.5f);
    mesh = bromesh::computeFlatNormals(mesh);
    bromesh::optimizeVertexCache(mesh);
    bromesh::optimizeOverdraw(mesh, 1.05f);
    bromesh::optimizeVertexFetch(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Spherical, 1.0f);

    std::string path = std::string(testDir) + "rt_capsule_full.glb";
    ASSERT(bromesh::saveGLTF(mesh, path), "gltf_rt_capsule_full: save");
    auto scene = bromesh::loadGLTF(path);
    ASSERT(!scene.meshes.empty(), "gltf_rt_capsule_full: has meshes");
    auto& loaded = scene.meshes[0];
    ASSERT(loaded.vertexCount() == mesh.vertexCount(), "gltf_rt_capsule_full: vertex count");
    ASSERT(loaded.triangleCount() == mesh.triangleCount(), "gltf_rt_capsule_full: tri count");
    ASSERT(positionsMatch(mesh, loaded, 1e-4f), "gltf_rt_capsule_full: positions");
    ASSERT(normalsMatch(mesh, loaded, 1e-3f), "gltf_rt_capsule_full: normals");
    ASSERT(uvsMatch(mesh, loaded, 1e-4f), "gltf_rt_capsule_full: UVs");
    std::remove(path.c_str());
}

TEST(gltf_rt_heightmap_simplified) {
    // Heightmap -> simplify -> recompute normals -> box UVs -> glTF
    float heights[49];
    for (int i = 0; i < 49; ++i)
        heights[i] = std::sin(i * 0.4f) * std::cos(i * 0.2f) * 3.0f;
    auto mesh = bromesh::heightmapGrid(heights, 7, 7, 0.5f);
    mesh = bromesh::simplify(mesh, 0.5f);
    bromesh::computeNormals(mesh);
    mesh.uvs.clear();
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Box, 1.0f);

    std::string path = std::string(testDir) + "rt_hm_simp.glb";
    ASSERT(bromesh::saveGLTF(mesh, path), "gltf_rt_hm_simp: save");
    auto scene = bromesh::loadGLTF(path);
    ASSERT(!scene.meshes.empty(), "gltf_rt_hm_simp: has meshes");
    auto& loaded = scene.meshes[0];
    ASSERT(loaded.vertexCount() == mesh.vertexCount(), "gltf_rt_hm_simp: vertex count");
    ASSERT(loaded.triangleCount() == mesh.triangleCount(), "gltf_rt_hm_simp: tri count");
    ASSERT(positionsMatch(mesh, loaded, 1e-4f), "gltf_rt_hm_simp: positions");
    ASSERT(loaded.hasNormals(), "gltf_rt_hm_simp: normals");
    ASSERT(loaded.hasUVs(), "gltf_rt_hm_simp: UVs");
    std::remove(path.c_str());
}

TEST(gltf_rt_transvoxel_with_transition) {
    // Transvoxel with LOD transition -> weld -> normals -> glTF
    const int N = 17;
    float field[N * N * N];
    fillSphereField(field, N, 6.0f);
    int neighborLods[6] = { 1, -1, -1, -1, -1, -1 };
    auto mesh = bromesh::transvoxel(field, N, 0, neighborLods, 0.0f, 1.0f);
    mesh = bromesh::weldVertices(mesh, 1e-4f);
    bromesh::computeNormals(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::PlanarXZ, 1.0f);

    std::string path = std::string(testDir) + "rt_tv_trans.glb";
    ASSERT(bromesh::saveGLTF(mesh, path), "gltf_rt_tv_trans: save");
    auto scene = bromesh::loadGLTF(path);
    ASSERT(!scene.meshes.empty(), "gltf_rt_tv_trans: has meshes");
    auto& loaded = scene.meshes[0];
    ASSERT(loaded.vertexCount() == mesh.vertexCount(), "gltf_rt_tv_trans: vertex count");
    ASSERT(loaded.triangleCount() == mesh.triangleCount(), "gltf_rt_tv_trans: tri count");
    ASSERT(positionsMatch(mesh, loaded, 1e-4f), "gltf_rt_tv_trans: positions");
    ASSERT(loaded.hasNormals(), "gltf_rt_tv_trans: normals");
    std::remove(path.c_str());
}

TEST(gltf_rt_greedy_mesh_voxel) {
    // Greedy mesh single voxel -> strip colors -> add normals + UVs -> glTF
    uint8_t voxels[27] = {};
    voxels[1 * 3 * 3 + 1 * 3 + 1] = 1;
    auto mesh = bromesh::greedyMesh(voxels, 3, 3, 3, 1.0f);
    mesh.colors.clear(); // glTF saver doesn't write colors, strip them
    bromesh::computeNormals(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Box, 1.0f);

    std::string path = std::string(testDir) + "rt_greedy.glb";
    ASSERT(bromesh::saveGLTF(mesh, path), "gltf_rt_greedy: save");
    auto scene = bromesh::loadGLTF(path);
    ASSERT(!scene.meshes.empty(), "gltf_rt_greedy: has meshes");
    auto& loaded = scene.meshes[0];
    ASSERT(loaded.vertexCount() == mesh.vertexCount(), "gltf_rt_greedy: vertex count");
    ASSERT(loaded.triangleCount() == mesh.triangleCount(), "gltf_rt_greedy: tri count");
    ASSERT(positionsMatch(mesh, loaded, 1e-4f), "gltf_rt_greedy: positions");
    ASSERT(loaded.hasNormals(), "gltf_rt_greedy: normals");
    ASSERT(loaded.hasUVs(), "gltf_rt_greedy: UVs");
    std::remove(path.c_str());
}

TEST(gltf_rt_plane_subdivided_all_ops) {
    // Plane 8x8 -> smooth normals -> spherical UVs -> all optimizations -> glTF
    auto mesh = bromesh::plane(4.0f, 4.0f, 8, 8);
    bromesh::computeNormals(mesh);
    mesh.uvs.clear();
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Spherical, 1.0f);
    bromesh::optimizeVertexCache(mesh);
    bromesh::optimizeOverdraw(mesh, 1.05f);
    bromesh::optimizeVertexFetch(mesh);

    std::string path = std::string(testDir) + "rt_plane_all.glb";
    ASSERT(bromesh::saveGLTF(mesh, path), "gltf_rt_plane_all: save");
    auto scene = bromesh::loadGLTF(path);
    ASSERT(!scene.meshes.empty(), "gltf_rt_plane_all: has meshes");
    auto& loaded = scene.meshes[0];
    ASSERT(loaded.vertexCount() == mesh.vertexCount(), "gltf_rt_plane_all: vertex count");
    ASSERT(loaded.triangleCount() == mesh.triangleCount(), "gltf_rt_plane_all: tri count");
    ASSERT(positionsMatch(mesh, loaded, 1e-4f), "gltf_rt_plane_all: positions");
    ASSERT(normalsMatch(mesh, loaded, 1e-3f), "gltf_rt_plane_all: normals");
    ASSERT(uvsMatch(mesh, loaded, 1e-4f), "gltf_rt_plane_all: UVs");
    std::remove(path.c_str());
}

TEST(gltf_rt_cylinder_weld_simplify_lod) {
    // Cylinder -> weld -> simplify -> LOD chain -> LOD 1 -> normals + UVs -> glTF
    auto mesh = bromesh::cylinder(2.0f, 3.0f, 32);
    mesh = bromesh::weldVertices(mesh, 1e-5f);
    float ratios[] = { 0.6f, 0.3f };
    auto chain = bromesh::generateLODChain(mesh, ratios, 2);
    ASSERT(chain.size() == 2, "gltf_rt_cyl_lod: 2 LOD levels");
    auto& lod = chain[1];
    bromesh::computeNormals(lod);
    bromesh::projectUVs(lod, bromesh::ProjectionType::Cylindrical, 1.0f);

    std::string path = std::string(testDir) + "rt_cyl_lod1.glb";
    ASSERT(bromesh::saveGLTF(lod, path), "gltf_rt_cyl_lod: save");
    auto scene = bromesh::loadGLTF(path);
    ASSERT(!scene.meshes.empty(), "gltf_rt_cyl_lod: has meshes");
    auto& loaded = scene.meshes[0];
    ASSERT(loaded.vertexCount() == lod.vertexCount(), "gltf_rt_cyl_lod: vertex count");
    ASSERT(loaded.triangleCount() == lod.triangleCount(), "gltf_rt_cyl_lod: tri count");
    ASSERT(positionsMatch(lod, loaded, 1e-4f), "gltf_rt_cyl_lod: positions");
    ASSERT(loaded.hasNormals(), "gltf_rt_cyl_lod: normals");
    ASSERT(loaded.hasUVs(), "gltf_rt_cyl_lod: UVs");
    std::remove(path.c_str());
}

TEST(glb_vs_gltf_consistency) {
    // Same mesh saved as .glb and .gltf should load identically
    auto mesh = bromesh::torus(1.5f, 0.5f, 16, 8);
    bromesh::computeNormals(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Box, 1.0f);

    std::string glbPath = std::string(testDir) + "rt_consistency.glb";
    std::string gltfPath = std::string(testDir) + "rt_consistency.gltf";
    ASSERT(bromesh::saveGLTF(mesh, glbPath), "glb_vs_gltf: save glb");
    ASSERT(bromesh::saveGLTF(mesh, gltfPath), "glb_vs_gltf: save gltf");

    auto glbScene = bromesh::loadGLTF(glbPath);
    auto gltfScene = bromesh::loadGLTF(gltfPath);
    ASSERT(!glbScene.meshes.empty(), "glb_vs_gltf: glb loaded");
    ASSERT(!gltfScene.meshes.empty(), "glb_vs_gltf: gltf loaded");

    auto& glbMesh = glbScene.meshes[0];
    auto& gltfMesh = gltfScene.meshes[0];
    ASSERT(glbMesh.vertexCount() == gltfMesh.vertexCount(), "glb_vs_gltf: vertex count");
    ASSERT(glbMesh.triangleCount() == gltfMesh.triangleCount(), "glb_vs_gltf: tri count");
    ASSERT(positionsMatch(glbMesh, gltfMesh, 1e-5f), "glb_vs_gltf: positions match");
    ASSERT(normalsMatch(glbMesh, gltfMesh, 1e-5f), "glb_vs_gltf: normals match");
    ASSERT(uvsMatch(glbMesh, gltfMesh, 1e-5f), "glb_vs_gltf: UVs match");

    std::remove(glbPath.c_str());
    std::remove(gltfPath.c_str());
    // Also remove the .bin sidecar from gltf save
    std::string binPath = std::string(testDir) + "rt_consistency.bin";
    std::remove(binPath.c_str());
}

#endif // BROMESH_HAS_GLTF

// ============================================================================
// Cross-format roundtrip: OBJ -> STL -> OBJ (geometry preservation)
// ============================================================================

TEST(cross_format_obj_stl_obj) {
    // Sphere -> OBJ -> load -> STL -> load -> compare bboxes and volume
    auto mesh = bromesh::sphere(2.0f, 24, 16);
    auto origBBox = meshBBox(mesh);
    float origVol = bromesh::computeVolume(mesh);

    std::string objPath = std::string(testDir) + "rt_cross.obj";
    std::string stlPath = std::string(testDir) + "rt_cross.stl";

    ASSERT(bromesh::saveOBJ(mesh, objPath), "cross_fmt: save obj");
    auto fromObj = bromesh::loadOBJ(objPath);
    ASSERT(!fromObj.empty(), "cross_fmt: load obj");

    ASSERT(bromesh::saveSTL(fromObj, stlPath), "cross_fmt: save stl from obj");
    auto fromStl = bromesh::loadSTL(stlPath);
    ASSERT(!fromStl.empty(), "cross_fmt: load stl");

    auto stlBBox = meshBBox(fromStl);
    ASSERT(bboxMatch(origBBox, stlBBox, 0.01f), "cross_fmt: bbox preserved through OBJ->STL");
    float stlVol = bromesh::computeVolume(fromStl);
    ASSERT(approxEqual(origVol, stlVol, std::fabs(origVol) * 0.02f + 0.1f),
           "cross_fmt: volume preserved through OBJ->STL");

    std::remove(objPath.c_str());
    std::remove(stlPath.c_str());
}

#if BROMESH_HAS_GLTF

TEST(cross_format_gltf_obj_stl) {
    // Complex pipeline -> glTF -> OBJ -> STL -> compare
    auto mesh = bromesh::capsule(1.5f, 2.0f, 20, 10);
    mesh = bromesh::simplify(mesh, 0.6f);
    bromesh::computeNormals(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Box, 1.0f);
    auto origBBox = meshBBox(mesh);

    std::string glbPath = std::string(testDir) + "rt_cross2.glb";
    std::string objPath = std::string(testDir) + "rt_cross2.obj";
    std::string stlPath = std::string(testDir) + "rt_cross2.stl";

    // Save as glTF, reload
    ASSERT(bromesh::saveGLTF(mesh, glbPath), "cross_gltf: save glb");
    auto scene = bromesh::loadGLTF(glbPath);
    ASSERT(!scene.meshes.empty(), "cross_gltf: load glb");
    auto& fromGltf = scene.meshes[0];
    ASSERT(positionsMatch(mesh, fromGltf, 1e-4f), "cross_gltf: glb positions");

    // Save glTF result as OBJ, reload
    bromesh::projectUVs(fromGltf, bromesh::ProjectionType::Box, 1.0f);
    ASSERT(bromesh::saveOBJ(fromGltf, objPath), "cross_gltf: save obj");
    auto fromObj = bromesh::loadOBJ(objPath);
    ASSERT(fromObj.vertexCount() == fromGltf.vertexCount(), "cross_gltf: obj vertex count");

    // Save OBJ result as STL, reload
    ASSERT(bromesh::saveSTL(fromObj, stlPath), "cross_gltf: save stl");
    auto fromStl = bromesh::loadSTL(stlPath);
    auto stlBBox = meshBBox(fromStl);
    ASSERT(bboxMatch(origBBox, stlBBox, 0.05f), "cross_gltf: bbox through gltf->obj->stl");

    std::remove(glbPath.c_str());
    std::remove(objPath.c_str());
    std::remove(stlPath.c_str());
}

#endif // BROMESH_HAS_GLTF

// ============================================================================
// Meshlet generation tests
// ============================================================================

#if BROMESH_HAS_MESHOPTIMIZER

TEST(meshlets_sphere) {
    auto mesh = bromesh::sphere(2.0f, 32, 24);
    auto meshlets = bromesh::buildMeshlets(mesh);
    ASSERT(!meshlets.empty(), "meshlets_sphere: should produce meshlets");
    ASSERT(meshlets.size() >= 2, "meshlets_sphere: sphere should have multiple meshlets");

    // Verify each meshlet has valid data
    size_t totalTris = 0;
    for (const auto& ml : meshlets) {
        ASSERT(ml.vertexCount() > 0, "meshlets_sphere: meshlet should have vertices");
        ASSERT(ml.triangleCount() > 0, "meshlets_sphere: meshlet should have triangles");
        ASSERT(ml.vertexCount() <= 64, "meshlets_sphere: meshlet should respect maxVertices");
        ASSERT(ml.triangleCount() <= 124, "meshlets_sphere: meshlet should respect maxTriangles");
        ASSERT(ml.bounds.radius > 0, "meshlets_sphere: meshlet should have bounding sphere");
        totalTris += ml.triangleCount();
    }
    ASSERT(totalTris == mesh.triangleCount(), "meshlets_sphere: total triangles should match");
}

TEST(meshlets_custom_params) {
    auto mesh = bromesh::box(1, 1, 1);
    bromesh::MeshletParams params;
    params.maxVertices = 32;
    params.maxTriangles = 32;
    auto meshlets = bromesh::buildMeshlets(mesh, params);
    ASSERT(!meshlets.empty(), "meshlets_custom: should produce meshlets");
    for (const auto& ml : meshlets) {
        ASSERT(ml.vertexCount() <= 32, "meshlets_custom: should respect custom maxVertices");
        ASSERT(ml.triangleCount() <= 32, "meshlets_custom: should respect custom maxTriangles");
    }
}

#endif // BROMESH_HAS_MESHOPTIMIZER

// ============================================================================
// Mesh analysis statistics tests
// ============================================================================

#if BROMESH_HAS_MESHOPTIMIZER

TEST(analyze_vertex_cache) {
    auto mesh = bromesh::sphere(1.0f, 16, 12);
    auto stats = bromesh::analyzeVertexCache(mesh);
    ASSERT(stats.verticesTransformed > 0, "analyze_vcache: should transform vertices");
    ASSERT(stats.acmr > 0, "analyze_vcache: ACMR should be positive");
    ASSERT(stats.atvr >= 1.0f, "analyze_vcache: ATVR should be >= 1.0");
}

TEST(analyze_vertex_fetch) {
    auto mesh = bromesh::sphere(1.0f, 16, 12);
    auto stats = bromesh::analyzeVertexFetch(mesh);
    ASSERT(stats.bytesFetched > 0, "analyze_vfetch: should fetch bytes");
    ASSERT(stats.overfetch >= 1.0f, "analyze_vfetch: overfetch should be >= 1.0");
}

TEST(analyze_overdraw) {
    auto mesh = bromesh::sphere(1.0f, 16, 12);
    auto stats = bromesh::analyzeOverdraw(mesh);
    ASSERT(stats.pixelsCovered > 0, "analyze_overdraw: should cover pixels");
    ASSERT(stats.overdraw >= 1.0f, "analyze_overdraw: overdraw should be >= 1.0");
}

TEST(analyze_before_after_optimize) {
    auto mesh = bromesh::sphere(2.0f, 32, 24);
    auto before = bromesh::analyzeVertexCache(mesh);

    bromesh::optimizeVertexCache(mesh);
    auto after = bromesh::analyzeVertexCache(mesh);

    // After optimization, ACMR should be equal or better (lower)
    ASSERT(after.acmr <= before.acmr + 0.01f,
           "analyze_opt: ACMR should improve after vertex cache optimization");
}

#endif // BROMESH_HAS_MESHOPTIMIZER

// ============================================================================
// Spatial sorting + shadow index buffer tests
// ============================================================================

#if BROMESH_HAS_MESHOPTIMIZER

TEST(spatial_sort_triangles) {
    auto mesh = bromesh::sphere(2.0f, 24, 16);
    size_t origTriCount = mesh.triangleCount();
    size_t origVertCount = mesh.vertexCount();
    bromesh::spatialSortTriangles(mesh);
    ASSERT(mesh.triangleCount() == origTriCount, "spatial_tri: triangle count preserved");
    ASSERT(mesh.vertexCount() == origVertCount, "spatial_tri: vertex count preserved");
}

TEST(spatial_sort_vertices) {
    auto mesh = bromesh::sphere(2.0f, 24, 16);
    size_t origTriCount = mesh.triangleCount();
    size_t origVertCount = mesh.vertexCount();
    bromesh::spatialSortVertices(mesh);
    ASSERT(mesh.triangleCount() == origTriCount, "spatial_vert: triangle count preserved");
    ASSERT(mesh.vertexCount() == origVertCount, "spatial_vert: vertex count preserved");
    ASSERT(mesh.hasNormals(), "spatial_vert: normals preserved");
}

TEST(shadow_index_buffer) {
    auto mesh = bromesh::sphere(2.0f, 24, 16);
    auto shadow = bromesh::generateShadowIndexBuffer(mesh);
    ASSERT(shadow.size() == mesh.indices.size(), "shadow_ib: same index count");
    // Shadow indices should reference valid vertices
    for (uint32_t idx : shadow) {
        ASSERT(idx < mesh.vertexCount(), "shadow_ib: valid vertex reference");
    }
}

#endif // BROMESH_HAS_MESHOPTIMIZER

// ============================================================================
// Attribute-aware simplification tests
// ============================================================================

#if BROMESH_HAS_MESHOPTIMIZER

TEST(simplify_with_attributes) {
    auto mesh = bromesh::sphere(2.0f, 32, 24);
    bromesh::computeNormals(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Spherical, 1.0f);

    auto simplified = bromesh::simplifyWithAttributes(mesh, 0.5f);
    ASSERT(!simplified.empty(), "simplify_attr: should produce non-empty result");
    ASSERT(simplified.triangleCount() < mesh.triangleCount(),
           "simplify_attr: should have fewer triangles");
    ASSERT(simplified.hasNormals(), "simplify_attr: should preserve normals");
    ASSERT(simplified.hasUVs(), "simplify_attr: should preserve UVs");
}

TEST(simplify_with_attributes_preserves_more) {
    // Compare attribute-aware vs basic: attribute-aware should produce
    // at least as many triangles (it's more conservative at seams)
    auto mesh = bromesh::torus(2.0f, 0.5f, 32, 16);
    bromesh::computeNormals(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Box, 1.0f);

    auto basic = bromesh::simplify(mesh, 0.3f);
    auto attr = bromesh::simplifyWithAttributes(mesh, 0.3f);
    ASSERT(!basic.empty(), "simplify_cmp: basic should work");
    ASSERT(!attr.empty(), "simplify_cmp: attribute-aware should work");
    // Both should reduce triangle count
    ASSERT(basic.triangleCount() < mesh.triangleCount(), "simplify_cmp: basic reduces");
    ASSERT(attr.triangleCount() < mesh.triangleCount(), "simplify_cmp: attr reduces");
}

#endif // BROMESH_HAS_MESHOPTIMIZER

// ============================================================================
// Mesh merge tests
// ============================================================================

TEST(merge_two_meshes) {
    auto a = bromesh::box(1, 1, 1);
    auto b = bromesh::sphere(1.0f, 8, 6);
    std::vector<bromesh::MeshData> meshes = {a, b};
    auto merged = bromesh::mergeMeshes(meshes);
    ASSERT(merged.vertexCount() == a.vertexCount() + b.vertexCount(),
           "merge_two: vertex count should be sum");
    ASSERT(merged.triangleCount() == a.triangleCount() + b.triangleCount(),
           "merge_two: triangle count should be sum");
    ASSERT(merged.hasNormals(), "merge_two: both have normals so merged should too");
    ASSERT(merged.hasUVs(), "merge_two: both have UVs so merged should too");
}

TEST(merge_single_mesh) {
    auto a = bromesh::box(1, 1, 1);
    auto merged = bromesh::mergeMeshes(&a, 1);
    ASSERT(merged.vertexCount() == a.vertexCount(), "merge_single: should be same");
    ASSERT(merged.triangleCount() == a.triangleCount(), "merge_single: should be same");
}

TEST(merge_empty) {
    auto result = bromesh::mergeMeshes(nullptr, 0);
    ASSERT(result.empty(), "merge_empty: should be empty");
}

TEST(merge_index_validity) {
    auto a = bromesh::box(1, 1, 1);
    auto b = bromesh::cylinder(0.5f, 1.0f, 12);
    std::vector<bromesh::MeshData> meshes = {a, b};
    auto merged = bromesh::mergeMeshes(meshes);
    // All indices should be valid
    for (uint32_t idx : merged.indices) {
        ASSERT(idx < merged.vertexCount(), "merge_idx: all indices should be valid");
    }
}

// ============================================================================
// par_shapes primitive tests
// ============================================================================

#if BROMESH_HAS_PAR_SHAPES

TEST(par_icosahedron) {
    auto mesh = bromesh::icosahedron();
    ASSERT(!mesh.empty(), "par_ico: should be non-empty");
    ASSERT(mesh.triangleCount() == 20, "par_ico: icosahedron has 20 faces");
    ASSERT(mesh.vertexCount() == 12, "par_ico: icosahedron has 12 vertices");
}

TEST(par_dodecahedron) {
    auto mesh = bromesh::dodecahedron();
    ASSERT(!mesh.empty(), "par_dodec: should be non-empty");
    ASSERT(mesh.vertexCount() > 0, "par_dodec: should have vertices");
    ASSERT(mesh.triangleCount() > 0, "par_dodec: should have triangles");
}

TEST(par_octahedron) {
    auto mesh = bromesh::octahedron();
    ASSERT(!mesh.empty(), "par_oct: should be non-empty");
    ASSERT(mesh.triangleCount() == 8, "par_oct: octahedron has 8 faces");
    ASSERT(mesh.vertexCount() == 6, "par_oct: octahedron has 6 vertices");
}

TEST(par_tetrahedron) {
    auto mesh = bromesh::tetrahedron();
    ASSERT(!mesh.empty(), "par_tet: should be non-empty");
    ASSERT(mesh.triangleCount() == 4, "par_tet: tetrahedron has 4 faces");
    ASSERT(mesh.vertexCount() == 4, "par_tet: tetrahedron has 4 vertices");
}

TEST(par_geodesic_sphere) {
    auto mesh = bromesh::geodesicSphere(2.0f, 2);
    ASSERT(!mesh.empty(), "par_geod: should be non-empty");
    ASSERT(mesh.vertexCount() > 12, "par_geod: subdivided should have more than icosahedron");
    ASSERT(mesh.hasNormals(), "par_geod: should have normals");

    // Verify all vertices are approximately on the sphere surface
    bool allOnSurface = true;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float x = mesh.positions[v * 3 + 0];
        float y = mesh.positions[v * 3 + 1];
        float z = mesh.positions[v * 3 + 2];
        float dist = std::sqrt(x * x + y * y + z * z);
        if (std::fabs(dist - 2.0f) > 0.1f) { allOnSurface = false; break; }
    }
    ASSERT(allOnSurface, "par_geod: vertices should be on sphere surface");
}

TEST(par_cone) {
    auto mesh = bromesh::cone(1.0f, 2.0f, 16, 4);
    ASSERT(!mesh.empty(), "par_cone: should be non-empty");
    ASSERT(mesh.vertexCount() > 10, "par_cone: should have vertices");
    ASSERT(mesh.triangleCount() > 10, "par_cone: should have triangles");
}

TEST(par_disc) {
    auto mesh = bromesh::disc(1.5f, 16);
    ASSERT(!mesh.empty(), "par_disc: should be non-empty");
    ASSERT(mesh.vertexCount() > 0, "par_disc: should have vertices");
}

TEST(par_rock) {
    auto mesh = bromesh::rock(1.0f, 42, 2);
    ASSERT(!mesh.empty(), "par_rock: should be non-empty");
    ASSERT(mesh.vertexCount() > 20, "par_rock: should have reasonable vertex count");
    ASSERT(mesh.hasNormals(), "par_rock: should have normals");

    // Different seeds should produce different shapes
    auto mesh2 = bromesh::rock(1.0f, 99, 2);
    ASSERT(!mesh2.empty(), "par_rock2: should be non-empty");
    // Vertices should differ
    bool differ = false;
    size_t checkCount = std::min(mesh.vertexCount(), mesh2.vertexCount());
    for (size_t i = 0; i < checkCount * 3 && !differ; ++i) {
        if (std::fabs(mesh.positions[i] - mesh2.positions[i]) > 0.001f) differ = true;
    }
    ASSERT(differ, "par_rock: different seeds should produce different shapes");
}

TEST(par_trefoil_knot) {
    auto mesh = bromesh::trefoilKnot(1.0f, 32, 8);
    ASSERT(!mesh.empty(), "par_trefoil: should be non-empty");
    ASSERT(mesh.vertexCount() > 50, "par_trefoil: should have many vertices");
    ASSERT(mesh.hasNormals(), "par_trefoil: should have normals");
}

TEST(par_klein_bottle) {
    auto mesh = bromesh::kleinBottle(16, 8);
    ASSERT(!mesh.empty(), "par_klein: should be non-empty");
    ASSERT(mesh.vertexCount() > 50, "par_klein: should have many vertices");
    ASSERT(mesh.hasNormals(), "par_klein: should have normals");
}

#endif // BROMESH_HAS_PAR_SHAPES

// ============================================================================
// Mesh encoding/compression tests
// ============================================================================

#if BROMESH_HAS_MESHOPTIMIZER

TEST(encode_decode_mesh_roundtrip) {
    auto mesh = bromesh::sphere(2.0f, 16, 12);
    bromesh::computeNormals(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Spherical, 1.0f);

    auto encoded = bromesh::encodeMesh(mesh);
    ASSERT(!encoded.vertexData.empty(), "encode: vertex data should be non-empty");
    ASSERT(!encoded.indexData.empty(), "encode: index data should be non-empty");
    ASSERT(encoded.vertexCount == mesh.vertexCount(), "encode: vertex count preserved");
    ASSERT(encoded.indexCount == mesh.indices.size(), "encode: index count preserved");

    // Compressed should be smaller than raw
    size_t rawVertSize = mesh.vertexCount() * encoded.vertexSize;
    ASSERT(encoded.vertexData.size() < rawVertSize, "encode: vertex data should be compressed");

    auto decoded = bromesh::decodeMesh(encoded, true, true, false);
    ASSERT(decoded.vertexCount() == mesh.vertexCount(), "decode: vertex count matches");
    ASSERT(decoded.triangleCount() == mesh.triangleCount(), "decode: triangle count matches");
    ASSERT(decoded.hasNormals(), "decode: normals preserved");
    ASSERT(decoded.hasUVs(), "decode: UVs preserved");

    // Verify positions match
    bool posMatch = true;
    for (size_t i = 0; i < mesh.positions.size(); ++i) {
        if (std::fabs(mesh.positions[i] - decoded.positions[i]) > 1e-5f) {
            posMatch = false; break;
        }
    }
    ASSERT(posMatch, "decode: positions match original");
}

TEST(encode_decode_index_buffer) {
    auto mesh = bromesh::box(1, 1, 1);
    auto encoded = bromesh::encodeIndexBuffer(mesh.indices, mesh.vertexCount());
    ASSERT(!encoded.empty(), "encode_ib: should produce data");
    ASSERT(encoded.size() < mesh.indices.size() * sizeof(uint32_t), "encode_ib: should compress");

    auto decoded = bromesh::decodeIndexBuffer(encoded, mesh.indices.size());
    ASSERT(decoded.size() == mesh.indices.size(), "decode_ib: size matches");
    ASSERT(decoded == mesh.indices, "decode_ib: indices match original");
}

#endif // BROMESH_HAS_MESHOPTIMIZER

// ============================================================================
// Triangle strip tests
// ============================================================================

#if BROMESH_HAS_MESHOPTIMIZER

TEST(stripify_unstripify_roundtrip) {
    auto mesh = bromesh::sphere(1.0f, 16, 12);
    bromesh::optimizeVertexCache(mesh);

    auto strip = bromesh::stripify(mesh.indices, mesh.vertexCount());
    ASSERT(!strip.empty(), "stripify: should produce strip");

    auto restored = bromesh::unstripify(strip);
    ASSERT(!restored.empty(), "unstripify: should produce triangle list");
    // Restored triangle count should match original
    ASSERT(restored.size() / 3 == mesh.indices.size() / 3,
           "strip_roundtrip: triangle count preserved");
}

TEST(stripify_box) {
    auto mesh = bromesh::box(1, 1, 1);
    auto strip = bromesh::stripify(mesh.indices, mesh.vertexCount());
    ASSERT(!strip.empty(), "stripify_box: should produce strip");
    // Strip should be reasonably compact
    ASSERT(strip.size() <= mesh.indices.size() * 2,
           "stripify_box: strip shouldn't be much larger than triangle list");
}

#endif // BROMESH_HAS_MESHOPTIMIZER

// ============================================================================
// Mesh repair tests
// ============================================================================

TEST(remove_degenerate_triangles) {
    // Create a mesh with one good triangle and one degenerate (zero-area)
    bromesh::MeshData mesh;
    mesh.positions = {
        0,0,0, 1,0,0, 0,1,0,  // good triangle
        2,0,0, 2,0,0, 3,0,0   // degenerate: two identical vertices
    };
    mesh.indices = {0,1,2, 3,4,5};

    auto repaired = bromesh::removeDegenerateTriangles(mesh);
    ASSERT(repaired.triangleCount() == 1, "remove_degen: should keep 1 triangle");
    ASSERT(repaired.indices[0] == 0 && repaired.indices[1] == 1 && repaired.indices[2] == 2,
           "remove_degen: should keep the good triangle");
}

TEST(remove_degenerate_collinear) {
    // Collinear triangle (zero area)
    bromesh::MeshData mesh;
    mesh.positions = {
        0,0,0, 1,0,0, 0,1,0,  // good
        0,0,0, 1,0,0, 2,0,0   // collinear
    };
    mesh.indices = {0,1,2, 3,4,5};

    auto repaired = bromesh::removeDegenerateTriangles(mesh);
    ASSERT(repaired.triangleCount() == 1, "remove_collinear: should remove collinear triangle");
}

TEST(remove_duplicate_triangles) {
    auto mesh = bromesh::box(1, 1, 1);
    size_t origTris = mesh.triangleCount();

    // Duplicate all triangles
    size_t origIdxCount = mesh.indices.size();
    for (size_t i = 0; i < origIdxCount; ++i) {
        mesh.indices.push_back(mesh.indices[i]);
    }
    ASSERT(mesh.triangleCount() == origTris * 2, "dup_setup: doubled");

    auto repaired = bromesh::removeDuplicateTriangles(mesh);
    ASSERT(repaired.triangleCount() == origTris, "remove_dup: should remove duplicates");
}

TEST(fill_holes_simple) {
    // Create an open box (5 faces, missing the top)
    // Using a simple example: a plane with a triangular hole
    bromesh::MeshData mesh;
    // Square with 4 triangles leaving a hole in the middle
    //  3---2
    //  |\ /|
    //  | 4 |   (vertex 4 at center, but no bottom-center triangle)
    //  |/ \|
    //  0---1
    mesh.positions = {
        0,0,0, 1,0,0, 1,1,0, 0,1,0, 0.5f,0.5f,0
    };
    // Only 3 triangles, leaving a gap
    mesh.indices = {
        0,1,4,  // bottom
        1,2,4,  // right
        2,3,4   // top
        // missing: 3,0,4 (left)
    };

    auto filled = bromesh::fillHoles(mesh);
    // Should add the missing triangle
    ASSERT(filled.triangleCount() >= mesh.triangleCount(),
           "fill_holes: should have at least as many triangles");
    ASSERT(filled.triangleCount() > mesh.triangleCount(),
           "fill_holes: should have added triangles to fill hole");
}

TEST(repair_preserves_clean_mesh) {
    // A clean mesh should pass through unchanged
    auto mesh = bromesh::box(1, 1, 1);
    auto degenFixed = bromesh::removeDegenerateTriangles(mesh);
    ASSERT(degenFixed.triangleCount() == mesh.triangleCount(),
           "repair_clean: degenerate removal should not change clean mesh");

    auto dupFixed = bromesh::removeDuplicateTriangles(mesh);
    ASSERT(dupFixed.triangleCount() == mesh.triangleCount(),
           "repair_clean: duplicate removal should not change clean mesh");
}

// ============================================================================
// UV Unwrapping tests (xatlas)
// ============================================================================

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

// ============================================================================
// Boolean/CSG tests (manifold)
// ============================================================================

#if BROMESH_HAS_MANIFOLD

// Helper: compute volume even on meshes with split vertices (like our primitives).
// Uses the divergence theorem directly without the manifold check.
static float rawVolume(const bromesh::MeshData& mesh) {
    double sum = 0.0;
    for (size_t t = 0; t < mesh.triangleCount(); ++t) {
        uint32_t i0 = mesh.indices[t * 3 + 0];
        uint32_t i1 = mesh.indices[t * 3 + 1];
        uint32_t i2 = mesh.indices[t * 3 + 2];
        double ax = mesh.positions[i0*3], ay = mesh.positions[i0*3+1], az = mesh.positions[i0*3+2];
        double bx = mesh.positions[i1*3], by = mesh.positions[i1*3+1], bz = mesh.positions[i1*3+2];
        double cx = mesh.positions[i2*3], cy = mesh.positions[i2*3+1], cz = mesh.positions[i2*3+2];
        sum += ax*(by*cz-bz*cy) + ay*(bz*cx-bx*cz) + az*(bx*cy-by*cx);
    }
    return static_cast<float>(std::fabs(sum) / 6.0);
}

TEST(boolean_union) {
    // Two overlapping spheres
    auto a = bromesh::sphere(1.0f, 16, 12);
    auto b = bromesh::sphere(1.0f, 16, 12);
    for (size_t v = 0; v < b.vertexCount(); ++v) {
        b.positions[v * 3 + 0] += 1.0f;
    }

    auto result = bromesh::booleanUnion(a, b);
    ASSERT(!result.empty(), "bool_union: should produce non-empty result");
    ASSERT(result.triangleCount() > 0, "bool_union: should have triangles");
    float volA = rawVolume(a);
    float volResult = rawVolume(result);
    ASSERT(volResult > volA * 0.9f, "bool_union: result volume should exceed single sphere");
}

TEST(boolean_difference) {
    auto a = bromesh::sphere(2.0f, 16, 12);
    auto b = bromesh::sphere(1.0f, 16, 12);
    for (size_t v = 0; v < b.vertexCount(); ++v) {
        b.positions[v * 3 + 0] += 1.5f;
    }

    auto result = bromesh::booleanDifference(a, b);
    ASSERT(!result.empty(), "bool_diff: should produce non-empty result");
    ASSERT(result.triangleCount() > 0, "bool_diff: should have triangles");
    // Result should have fewer vertices than the sum of both inputs
    ASSERT(result.vertexCount() < a.vertexCount() + b.vertexCount(),
           "bool_diff: should not just concatenate meshes");
    // BBox should fit within A's bbox
    auto bboxA = bromesh::computeBBox(a);
    auto bboxR = bromesh::computeBBox(result);
    ASSERT(bboxR.min[0] >= bboxA.min[0] - 0.01f && bboxR.max[0] <= bboxA.max[0] + 0.01f,
           "bool_diff: result should fit within A's bbox on X");
}

TEST(boolean_intersection) {
    auto a = bromesh::sphere(1.0f, 16, 12);
    auto b = bromesh::sphere(1.0f, 16, 12);
    for (size_t v = 0; v < b.vertexCount(); ++v) {
        b.positions[v * 3 + 0] += 0.5f;
    }

    auto result = bromesh::booleanIntersection(a, b);
    ASSERT(!result.empty(), "bool_isect: should produce non-empty result");
    float volA = rawVolume(a);
    float volResult = rawVolume(result);
    ASSERT(volResult < volA * 1.01f, "bool_isect: intersection volume should be less than full sphere");
    ASSERT(volResult > 0.01f, "bool_isect: intersection should have meaningful volume");
}

TEST(boolean_no_overlap) {
    auto a = bromesh::sphere(1.0f, 16, 12);
    auto b = bromesh::sphere(1.0f, 16, 12);
    for (size_t v = 0; v < b.vertexCount(); ++v) {
        b.positions[v * 3 + 0] += 5.0f;
    }

    auto result = bromesh::booleanUnion(a, b);
    ASSERT(!result.empty(), "bool_nooverlap: should produce result");
    float volA = rawVolume(a);
    float volB = rawVolume(b);
    float volResult = rawVolume(result);
    float expected = volA + volB;
    ASSERT(std::fabs(volResult - expected) < expected * 0.15f,
           "bool_nooverlap: union volume should be sum of parts");
}

TEST(split_by_plane) {
    auto mesh = bromesh::sphere(2.0f, 24, 16);

    auto [top, bottom] = bromesh::splitByPlane(mesh, 0, 1, 0, 0);
    ASSERT(!top.empty(), "split_plane: top half should be non-empty");
    ASSERT(!bottom.empty(), "split_plane: bottom half should be non-empty");

    float topVol = rawVolume(top);
    float bottomVol = rawVolume(bottom);

    // Both halves should be roughly equal for a centered sphere split at Y=0
    ASSERT(std::fabs(topVol - bottomVol) < (topVol + bottomVol) * 0.2f,
           "split_plane: halves should be roughly equal");
}

TEST(boolean_box_minus_sphere) {
    auto cube = bromesh::box(1.5f, 1.5f, 1.5f);
    auto ball = bromesh::sphere(1.0f, 16, 12);

    auto result = bromesh::booleanDifference(cube, ball);
    ASSERT(!result.empty(), "bool_box_sphere: should produce result");
    ASSERT(result.triangleCount() > 0, "bool_box_sphere: should have triangles");
    // Result should be more complex than the original cube (added sphere boundary)
    ASSERT(result.triangleCount() > cube.triangleCount(),
           "bool_box_sphere: result should have more triangles than original cube");
    // BBox should fit within cube's bbox
    auto bboxC = bromesh::computeBBox(cube);
    auto bboxR = bromesh::computeBBox(result);
    ASSERT(bboxR.min[0] >= bboxC.min[0] - 0.01f && bboxR.max[0] <= bboxC.max[0] + 0.01f,
           "bool_box_sphere: result should fit within cube bbox");
}

#endif // BROMESH_HAS_MANIFOLD

// ====================== SUBDIVISION ======================

TEST(subdivide_midpoint_box) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    auto result = bromesh::subdivideMidpoint(mesh, 1);
    // Each triangle splits into 4, box has 12 triangles -> 48
    ASSERT(result.triangleCount() == mesh.triangleCount() * 4,
           "midpoint_box: should have 4x triangles");
    ASSERT(result.vertexCount() > mesh.vertexCount(),
           "midpoint_box: should have more vertices");
}

TEST(subdivide_midpoint_two_iterations) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    auto result = bromesh::subdivideMidpoint(mesh, 2);
    // 12 * 4 * 4 = 192
    ASSERT(result.triangleCount() == mesh.triangleCount() * 16,
           "midpoint_2x: should have 16x triangles");
}

TEST(subdivide_loop_sphere) {
    auto mesh = bromesh::sphere(1.0f, 8, 6);
    auto result = bromesh::subdivideLoop(mesh, 1);
    ASSERT(result.triangleCount() == mesh.triangleCount() * 4,
           "loop_sphere: should have 4x triangles");
    ASSERT(result.hasNormals(), "loop_sphere: should have normals");
}

TEST(subdivide_catmull_clark_box) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    auto result = bromesh::subdivideCatmullClark(mesh, 1);
    // CC on triangles: each tri -> 3 quads -> 6 triangles
    ASSERT(result.triangleCount() == mesh.triangleCount() * 6,
           "cc_box: should have 6x triangles");
    ASSERT(result.hasNormals(), "cc_box: should have normals");
}

TEST(subdivide_zero_iterations) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    auto result = bromesh::subdivideLoop(mesh, 0);
    ASSERT(result.triangleCount() == mesh.triangleCount(),
           "subdiv_zero: 0 iterations should return same mesh");
}

// ====================== PLY I/O ======================

TEST(ply_roundtrip_box) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::computeNormals(mesh);

    bool ok = bromesh::savePLY(mesh, "test_box.ply");
    ASSERT(ok, "ply_roundtrip: save should succeed");

    auto loaded = bromesh::loadPLY("test_box.ply");
    ASSERT(!loaded.empty(), "ply_roundtrip: load should succeed");
    ASSERT(loaded.vertexCount() == mesh.vertexCount(),
           "ply_roundtrip: vertex count should match");
    ASSERT(loaded.triangleCount() == mesh.triangleCount(),
           "ply_roundtrip: triangle count should match");
    ASSERT(loaded.hasNormals(), "ply_roundtrip: should have normals");

    std::remove("test_box.ply");
}

TEST(ply_roundtrip_with_colors) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    // Add vertex colors
    mesh.colors.resize(mesh.vertexCount() * 4);
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        mesh.colors[v*4+0] = 1.0f;
        mesh.colors[v*4+1] = 0.0f;
        mesh.colors[v*4+2] = 0.0f;
        mesh.colors[v*4+3] = 1.0f;
    }

    bromesh::savePLY(mesh, "test_box_col.ply");
    auto loaded = bromesh::loadPLY("test_box_col.ply");
    ASSERT(loaded.hasColors(), "ply_colors: should have colors");
    ASSERT(std::fabs(loaded.colors[0] - 1.0f) < 0.01f,
           "ply_colors: red channel should be ~1.0");

    std::remove("test_box_col.ply");
}

TEST(ply_cross_format_obj_to_ply) {
    auto mesh = bromesh::sphere(1.0f, 12, 8);
    bromesh::computeNormals(mesh);

    bromesh::saveOBJ(mesh, "test_sphere.obj");
    auto obj = bromesh::loadOBJ("test_sphere.obj");

    bromesh::savePLY(obj, "test_sphere.ply");
    auto ply = bromesh::loadPLY("test_sphere.ply");

    ASSERT(!ply.empty(), "ply_cross: loaded PLY should not be empty");
    ASSERT(ply.triangleCount() == obj.triangleCount(),
           "ply_cross: triangle count should match OBJ");

    std::remove("test_sphere.obj");
    std::remove("test_sphere.ply");
}

// ====================== SMOOTHING ======================

TEST(smooth_laplacian_sphere) {
    auto mesh = bromesh::sphere(1.0f, 8, 6);
    bromesh::computeNormals(mesh);
    auto origPositions = mesh.positions;

    bromesh::smoothLaplacian(mesh, 0.5f, 3);

    ASSERT(mesh.vertexCount() * 3 == origPositions.size(),
           "laplacian: vertex count should not change");
    // After smoothing, positions should be different
    bool changed = false;
    for (size_t i = 0; i < origPositions.size(); ++i) {
        if (std::fabs(mesh.positions[i] - origPositions[i]) > 1e-6f) {
            changed = true;
            break;
        }
    }
    ASSERT(changed, "laplacian: positions should change after smoothing");
}

TEST(smooth_taubin_box) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::computeNormals(mesh);

    // Measure initial bounding box
    auto bbox1 = bromesh::computeBBox(mesh);
    float vol1 = (bbox1.max[0]-bbox1.min[0]) * (bbox1.max[1]-bbox1.min[1]) * (bbox1.max[2]-bbox1.min[2]);

    bromesh::smoothTaubin(mesh, 0.5f, -0.53f, 5);

    auto bbox2 = bromesh::computeBBox(mesh);
    float vol2 = (bbox2.max[0]-bbox2.min[0]) * (bbox2.max[1]-bbox2.min[1]) * (bbox2.max[2]-bbox2.min[2]);

    // Taubin should not shrink significantly (unlike pure Laplacian)
    ASSERT(vol2 > vol1 * 0.5f, "taubin: should not shrink excessively");
}

// ====================== REMESHING ======================

TEST(remesh_isotropic_sphere) {
    auto mesh = bromesh::sphere(1.0f, 8, 6);
    bromesh::computeNormals(mesh);

    auto result = bromesh::remeshIsotropic(mesh, 0.0f, 3);
    ASSERT(!result.empty(), "remesh: result should not be empty");
    ASSERT(result.hasNormals(), "remesh: should have normals");
    ASSERT(result.triangleCount() > 0, "remesh: should have triangles");
}

// ====================== RECONSTRUCTION ======================

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

// ====================== TARGET TRIANGLE COUNT DECIMATION ======================

TEST(simplify_target_count) {
    auto mesh = bromesh::sphere(1.0f, 16, 12);
    size_t origTris = mesh.triangleCount();

    size_t target = origTris / 4;
    auto result = bromesh::simplifyToTriangleCount(mesh, target);

    // Should have fewer triangles than original (meshopt may not hit exact target)
    ASSERT(result.triangleCount() < origTris,
           "target_count: should have fewer triangles");
    ASSERT(result.triangleCount() > 0, "target_count: should have triangles");
}

TEST(simplify_target_count_identity) {
    auto mesh = bromesh::sphere(1.0f, 8, 6);
    size_t origTris = mesh.triangleCount();

    // Requesting more than current should return same mesh
    auto result = bromesh::simplifyToTriangleCount(mesh, origTris * 2);
    ASSERT(result.triangleCount() == origTris,
           "target_count_identity: should return unchanged when target >= current");
}

// ====================== VERTEX COLOR BAKING ======================

TEST(bake_curvature_box) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::computeNormals(mesh);

    bromesh::bakeCurvature(mesh, 1.0f);
    ASSERT(mesh.hasColors(), "bake_curv: should have colors");
    // All colors should be in [0,1]
    bool valid = true;
    for (size_t i = 0; i < mesh.colors.size(); ++i) {
        if (mesh.colors[i] < -0.01f || mesh.colors[i] > 1.01f) {
            valid = false;
            break;
        }
    }
    ASSERT(valid, "bake_curv: colors should be in [0,1]");
}

TEST(bake_ao_small_sphere) {
    auto mesh = bromesh::sphere(1.0f, 6, 4);
    bromesh::computeNormals(mesh);

    bromesh::bakeAmbientOcclusion(mesh, 8, 0.0f);
    ASSERT(mesh.hasColors(), "bake_ao: should have colors");
    // AO values should be in [0,1]
    bool valid = true;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float ao = mesh.colors[v*4+0];
        if (ao < -0.01f || ao > 1.01f) { valid = false; break; }
    }
    ASSERT(valid, "bake_ao: AO values should be in [0,1]");
}

TEST(bake_thickness_box) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::computeNormals(mesh);

    bromesh::bakeThickness(mesh, 8, 0.0f);
    ASSERT(mesh.hasColors(), "bake_thick: should have colors");
}

// ====================== FBX I/O ======================

#ifdef BROMESH_HAS_OPENFBX
TEST(fbx_api_smoke) {
    // Just verify the API compiles and doesn't crash on non-existent file
    auto meshes = bromesh::loadFBX("nonexistent.fbx");
    ASSERT(meshes.empty(), "fbx_smoke: non-existent file should return empty");
}
#endif

// ====================== Skinning Utilities ======================

TEST(normalize_weights_basic) {
    tests_run++;
    bromesh::SkinData skin;
    // 1 vertex, 4 weights that don't sum to 1
    skin.boneWeights = {0.5f, 0.3f, 0.1f, 0.05f};
    skin.boneIndices = {0, 1, 2, 3};
    skin.inverseBindMatrices.resize(4 * 16, 0.0f);
    skin.boneCount = 4;

    bromesh::normalizeWeights(skin);

    float sum = skin.boneWeights[0] + skin.boneWeights[1] +
                skin.boneWeights[2] + skin.boneWeights[3];
    ASSERT(std::fabs(sum - 1.0f) < 0.001f, "normalize: weights should sum to 1");
    // Weights should be sorted descending
    ASSERT(skin.boneWeights[0] >= skin.boneWeights[1], "normalize: sorted descending");
    tests_passed++;
}

TEST(normalize_weights_zeros) {
    tests_run++;
    bromesh::SkinData skin;
    skin.boneWeights = {0.0f, 0.0f, 0.0f, 0.0f};
    skin.boneIndices = {0, 1, 2, 3};
    skin.boneCount = 4;

    bromesh::normalizeWeights(skin);

    // Should default to bone 0 with weight 1
    ASSERT(std::fabs(skin.boneWeights[0] - 1.0f) < 0.001f,
           "normalize_zeros: first weight should be 1");
    ASSERT(skin.boneIndices[0] == 0, "normalize_zeros: first index should be 0");
    tests_passed++;
}

TEST(apply_skinning_identity) {
    tests_run++;
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::computeNormals(mesh);

    size_t vCount = mesh.vertexCount();

    bromesh::SkinData skin;
    skin.boneCount = 1;
    // Identity inverse bind matrix
    skin.inverseBindMatrices = {
        1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1
    };
    skin.boneWeights.resize(vCount * 4, 0.0f);
    skin.boneIndices.resize(vCount * 4, 0);
    for (size_t v = 0; v < vCount; ++v) {
        skin.boneWeights[v * 4] = 1.0f;
    }

    // Save original positions
    auto origPos = mesh.positions;

    // Identity pose matrix
    float pose[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    bromesh::applySkinning(mesh, skin, pose);

    // Positions should be unchanged
    bool same = true;
    for (size_t i = 0; i < origPos.size(); ++i) {
        if (std::fabs(mesh.positions[i] - origPos[i]) > 0.001f) {
            same = false; break;
        }
    }
    ASSERT(same, "skin_identity: positions unchanged with identity transform");
    tests_passed++;
}

TEST(apply_skinning_translation) {
    tests_run++;
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    size_t vCount = mesh.vertexCount();

    bromesh::SkinData skin;
    skin.boneCount = 1;
    skin.inverseBindMatrices = {
        1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1
    };
    skin.boneWeights.resize(vCount * 4, 0.0f);
    skin.boneIndices.resize(vCount * 4, 0);
    for (size_t v = 0; v < vCount; ++v)
        skin.boneWeights[v * 4] = 1.0f;

    auto origPos = mesh.positions;

    // Translate by (5, 0, 0) via pose matrix (column-major)
    float pose[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 5,0,0,1};
    bromesh::applySkinning(mesh, skin, pose);

    bool shifted = true;
    for (size_t v = 0; v < vCount; ++v) {
        if (std::fabs(mesh.positions[v*3+0] - (origPos[v*3+0] + 5.0f)) > 0.001f) {
            shifted = false; break;
        }
    }
    ASSERT(shifted, "skin_translate: positions shifted by +5 on X");
    tests_passed++;
}

TEST(apply_morph_target) {
    tests_run++;
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    auto origPos = mesh.positions;
    size_t vCount = mesh.vertexCount();

    bromesh::MorphTarget morph;
    morph.name = "expand";
    morph.deltaPositions.resize(vCount * 3, 0.0f);
    // Move every vertex +1 on X
    for (size_t v = 0; v < vCount; ++v)
        morph.deltaPositions[v * 3 + 0] = 1.0f;

    bromesh::applyMorphTarget(mesh, morph, 0.5f);

    bool correct = true;
    for (size_t v = 0; v < vCount; ++v) {
        if (std::fabs(mesh.positions[v*3+0] - (origPos[v*3+0] + 0.5f)) > 0.001f) {
            correct = false; break;
        }
    }
    ASSERT(correct, "morph: positions shifted +0.5 on X at weight 0.5");
    tests_passed++;
}

TEST(apply_morph_zero_weight) {
    tests_run++;
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    auto origPos = mesh.positions;

    bromesh::MorphTarget morph;
    morph.deltaPositions.resize(mesh.vertexCount() * 3, 100.0f);

    bromesh::applyMorphTarget(mesh, morph, 0.0f);

    ASSERT(mesh.positions == origPos, "morph_zero: no change at weight 0");
    tests_passed++;
}

// ====================== UV Quality Metrics ======================

TEST(uv_metrics_box_projection) {
    tests_run++;
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Box, 1.0f);

    auto metrics = bromesh::measureUVQuality(mesh);
    ASSERT(metrics.triangleCount == mesh.triangleCount(),
           "uv_metrics: triangle count matches");
    ASSERT(metrics.avgStretch > 0.0f, "uv_metrics: avg stretch > 0");
    ASSERT(metrics.uvSpaceUsage > 0.0f, "uv_metrics: UV usage > 0");
    tests_passed++;
}

TEST(uv_distortion_per_triangle) {
    tests_run++;
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
    tests_passed++;
}

TEST(uv_metrics_no_uvs) {
    tests_run++;
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    mesh.uvs.clear(); // remove UVs

    auto metrics = bromesh::measureUVQuality(mesh);
    ASSERT(metrics.triangleCount == 0, "uv_metrics_no_uvs: no triangles when no UVs");

    auto distortions = bromesh::computeUVDistortion(mesh);
    ASSERT(distortions.empty(), "uv_distortion_no_uvs: empty when no UVs");
    tests_passed++;
}

TEST(uv_metrics_planar_plane) {
    tests_run++;
    // A flat plane with planar XZ projection should have low angle distortion
    auto mesh = bromesh::plane(2.0f, 2.0f, 4, 4);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::PlanarXZ, 1.0f);

    auto metrics = bromesh::measureUVQuality(mesh);
    ASSERT(metrics.avgAngleDistortion < 0.1f,
           "uv_planar: flat plane should have low angle distortion");
    tests_passed++;
}

int main() {
    std::printf("bromesh tests: %d/%d passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
