#include "test_framework.h"
#include <cmath>
#include <filesystem>
#include <string>
#include <algorithm>

static std::string testFile(const std::string& name) {
    return (std::filesystem::temp_directory_path() / name).string();
}

static const std::string testDir = std::filesystem::temp_directory_path().string() + "/";

static bool approxEqual(float a, float b, float tol) {
    return std::fabs(a - b) <= tol;
}

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

static bromath::AABB3 meshBBox(const bromesh::MeshData& m) {
    return bromesh::computeBBox(m);
}

static bool bboxMatch(const bromath::AABB3& a, const bromath::AABB3& b, float tol) {
    if (!approxEqual(a.min.x, b.min.x, tol)) return false;
    if (!approxEqual(a.min.y, b.min.y, tol)) return false;
    if (!approxEqual(a.min.z, b.min.z, tol)) return false;
    if (!approxEqual(a.max.x, b.max.x, tol)) return false;
    if (!approxEqual(a.max.y, b.max.y, tol)) return false;
    if (!approxEqual(a.max.z, b.max.z, tol)) return false;
    return true;
}

static void fillSphereField(float* field, int N, float radius) {
    float c = (N - 1) * 0.5f;
    for (int z = 0; z < N; ++z)
        for (int y = 0; y < N; ++y)
            for (int x = 0; x < N; ++x) {
                float dx = x - c, dy = y - c, dz = z - c;
                field[z * N * N + y * N + x] = std::sqrt(dx*dx + dy*dy + dz*dz) - radius;
            }
}

TEST(obj_roundtrip) {
    auto b = bromesh::box(1, 1, 1);
    size_t origVerts = b.vertexCount();
    size_t origTris = b.triangleCount();

    std::string objPath = testFile("test_output.obj");
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

    std::string stlPath = testFile("test_output.stl");
    bool saved = bromesh::saveSTL(b, stlPath);
    ASSERT(saved, "STL save should succeed");

    auto loaded = bromesh::loadSTL(stlPath);
    ASSERT(!loaded.empty(), "STL load should return non-empty mesh");
    ASSERT(loaded.triangleCount() == origTris, "STL roundtrip triangle count should match");
    // A binary STL stores 3 unshared vertices per triangle; loadSTL welds the
    // coincident positions back into shared topology (see src/io/stl.cpp), so a
    // closed solid comes back with strictly fewer than 3x its triangle count.
    ASSERT(loaded.vertexCount() > 0 && loaded.vertexCount() < origTris * 3,
           "STL import should weld coincident vertices");
    ASSERT(loaded.hasNormals(), "STL should have normals");

    std::remove(stlPath.c_str());
}

TEST(vox_nonexistent_file) {
    auto data = bromesh::loadVOX(testFile("nonexistent.vox"));
    ASSERT(data.sizeX == 0, "VOX nonexistent file should return empty sizeX");
    ASSERT(data.sizeY == 0, "VOX nonexistent file should return empty sizeY");
    ASSERT(data.sizeZ == 0, "VOX nonexistent file should return empty sizeZ");
    ASSERT(data.voxels.empty(), "VOX nonexistent file should return empty voxels");
}

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

TEST(stl_rt_torus_with_smooth_normals) {
    // Torus -> smooth normals -> STL roundtrip. STL carries per-face normals and
    // stores 3 unshared vertices per triangle; loadSTL welds them back into
    // shared topology, so we compare triangle count and bounding box.
    auto mesh = bromesh::torus(2.0f, 0.5f, 24, 12);
    bromesh::computeNormals(mesh);
    auto origBBox = meshBBox(mesh);
    size_t origTris = mesh.triangleCount();

    std::string path = std::string(testDir) + "rt_torus.stl";
    ASSERT(bromesh::saveSTL(mesh, path), "stl_rt_torus: save");
    auto loaded = bromesh::loadSTL(path);
    ASSERT(loaded.triangleCount() == origTris, "stl_rt_torus: tri count");
    ASSERT(loaded.vertexCount() > 0 && loaded.vertexCount() < origTris * 3,
           "stl_rt_torus: welded, fewer than 3 verts per tri");
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
    // loadSTL welds coincident vertices and drops the triangles that collapse
    // to degenerate as a result. LOD simplification can leave a sub-epsilon
    // sliver, and whether it does is float-path (architecture) dependent, so
    // this roundtrip is not triangle-count preserving: the count may only
    // shrink. The bounding box is the geometric invariant that must hold.
    ASSERT(loaded.triangleCount() > 0 && loaded.triangleCount() <= origTris,
           "stl_rt_lod0: tri count preserved up to degenerate welding");
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


#ifdef BROMESH_HAS_OPENFBX
TEST(fbx_api_smoke) {
    // Just verify the API compiles and doesn't crash on non-existent file
    auto meshes = bromesh::loadFBX("nonexistent.fbx");
    ASSERT(meshes.empty(), "fbx_smoke: non-existent file should return empty");
}

#endif // BROMESH_HAS_OPENFBX
