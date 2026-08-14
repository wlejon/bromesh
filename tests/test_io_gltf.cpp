#include "test_framework.h"
#include <cmath>
#include <filesystem>
#include <string>
#include <algorithm>
#include <cstring>

#if BROMESH_HAS_GLTF

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

TEST(gltf_animation_roundtrip_ground_truth) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::computeNormals(mesh);
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Box, 1.0f);

    bromesh::Skeleton skeleton;
    bromesh::Bone b0;
    b0.name = "root";
    b0.parent = -1;
    b0.localT[0] = 0.0f; b0.localT[1] = 0.0f; b0.localT[2] = 0.0f;
    b0.localR[0] = 0.0f; b0.localR[1] = 0.0f; b0.localR[2] = 0.0f; b0.localR[3] = 1.0f;
    b0.localS[0] = 1.0f; b0.localS[1] = 1.0f; b0.localS[2] = 1.0f;

    bromesh::Bone b1;
    b1.name = "bone1";
    b1.parent = 0;
    b1.localT[0] = 0.0f; b1.localT[1] = 1.0f; b1.localT[2] = 0.0f;
    b1.localR[0] = 0.0f; b1.localR[1] = 0.0f; b1.localR[2] = 0.0f; b1.localR[3] = 1.0f;
    b1.localS[0] = 1.0f; b1.localS[1] = 1.0f; b1.localS[2] = 1.0f;

    skeleton.bones.push_back(b0);
    skeleton.bones.push_back(b1);

    bromesh::SkinData skin;
    skin.boneCount = 2;
    skin.boneIndices.assign(mesh.vertexCount() * 4, 0u);
    skin.boneWeights.assign(mesh.vertexCount() * 4, 0.0f);
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        skin.boneWeights[v * 4 + 0] = 1.0f; // All weighted to bone 0
    }
    skin.inverseBindMatrices.assign(2 * 16, 0.0f);
    float idMat[16] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
    std::memcpy(&skin.inverseBindMatrices[0], idMat, sizeof(idMat));
    std::memcpy(&skin.inverseBindMatrices[16], idMat, sizeof(idMat));

    bromesh::Animation anim;
    anim.name = "TestAnim";
    anim.duration = 1.0f;

    // Channel 0: Translation on bone 0 (Linear)
    bromesh::AnimChannel ch0;
    ch0.boneIndex = 0;
    ch0.path = bromesh::AnimChannel::Path::Translation;
    ch0.interp = bromesh::AnimChannel::Interp::Linear;
    ch0.times = { 0.0f, 0.25f, 0.5f, 1.0f };
    ch0.values = {
        0.0f, 0.0f, 0.0f,
        1.0f, 0.0f, 0.0f,
        2.0f, 1.0f, 0.0f,
        0.0f, 2.0f, 1.0f
    };
    anim.channels.push_back(ch0);

    // Channel 1: Rotation on bone 1 (Linear)
    bromesh::AnimChannel ch1;
    ch1.boneIndex = 1;
    ch1.path = bromesh::AnimChannel::Path::Rotation;
    ch1.interp = bromesh::AnimChannel::Interp::Linear;
    ch1.times = { 0.0f, 0.5f, 1.0f };
    ch1.values = {
        0.0f, 0.0f, 0.0f, 1.0f,
        0.0f, 0.7071068f, 0.0f, 0.7071068f,
        0.0f, 1.0f, 0.0f, 0.0f
    };
    anim.channels.push_back(ch1);

    // Channel 2: Scale on bone 1 (Step)
    bromesh::AnimChannel ch2;
    ch2.boneIndex = 1;
    ch2.path = bromesh::AnimChannel::Path::Scale;
    ch2.interp = bromesh::AnimChannel::Interp::Step;
    ch2.times = { 0.0f, 0.5f, 1.0f };
    ch2.values = {
        1.0f, 1.0f, 1.0f,
        2.0f, 2.0f, 2.0f,
        1.5f, 1.5f, 1.5f
    };
    anim.channels.push_back(ch2);

    std::string path = std::string(testDir) + "rt_anim_ground_truth.glb";
    ASSERT(bromesh::saveGLTF(mesh, &skin, &skeleton, {anim}, path), "saveGLTF with animation");

    auto scene = bromesh::loadGLTF(path);
    ASSERT(scene.animations.size() == 1, "loaded scene has 1 animation");
    const auto& loadedAnim = scene.animations[0];
    ASSERT(loadedAnim.channels.size() == anim.channels.size(), "channel count matches");

    for (size_t c = 0; c < anim.channels.size(); ++c) {
        const auto& origCh = anim.channels[c];
        bool found = false;
        for (const auto& lCh : loadedAnim.channels) {
            if (lCh.boneIndex == origCh.boneIndex && lCh.path == origCh.path) {
                found = true;
                ASSERT(lCh.interp == origCh.interp, "interpolation mode matches");
                ASSERT(lCh.times.size() == origCh.times.size(), "timestamp count matches");
                for (size_t t = 0; t < origCh.times.size(); ++t) {
                    ASSERT(std::fabs(lCh.times[t] - origCh.times[t]) < 1e-4f, "timestamp matches to 1e-4");
                }
                ASSERT(lCh.values.size() == origCh.values.size(), "value count matches");
                for (size_t v = 0; v < origCh.values.size(); ++v) {
                    ASSERT(std::fabs(lCh.values[v] - origCh.values[v]) < 1e-4f, "keyframe value matches to 1e-4");
                }
                break;
            }
        }
        ASSERT(found, "matching channel found in loaded animation");
    }

    std::remove(path.c_str());
}

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
        ASSERT(positionsMatch(mesh, loaded, 1e-4f), "gltf_rt_box_uvs: positions");
        ASSERT(normalsMatch(mesh, loaded, 1e-4f), "gltf_rt_box_uvs: normals");
        ASSERT(uvsMatch(mesh, loaded, 1e-4f), "gltf_rt_box_uvs: UVs");
        std::remove(path.c_str());
    }
}

TEST(gltf_rt_sphere_simplified_cylindrical) {
    // Sphere -> simplify -> recompute normals -> cylindrical UVs -> glTF roundtrip
    auto mesh = bromesh::sphere(2.5f, 32, 24);
    mesh = bromesh::simplify(mesh, 0.4f);
    ASSERT(!mesh.empty(), "gltf_rt_sphere: simplified not empty");
    bromesh::computeNormals(mesh);
    mesh.uvs.clear();
    bromesh::projectUVs(mesh, bromesh::ProjectionType::Cylindrical, 1.0f);

    std::string path = std::string(testDir) + "rt_sphere_simp.glb";
    ASSERT(bromesh::saveGLTF(mesh, path), "gltf_rt_sphere: save");
    auto scene = bromesh::loadGLTF(path);
    ASSERT(!scene.meshes.empty(), "gltf_rt_sphere: has meshes");
    auto& loaded = scene.meshes[0];
    ASSERT(loaded.vertexCount() == mesh.vertexCount(), "gltf_rt_sphere: vertex count");
    ASSERT(loaded.triangleCount() == mesh.triangleCount(), "gltf_rt_sphere: tri count");
    ASSERT(positionsMatch(mesh, loaded, 1e-4f), "gltf_rt_sphere: positions");
    ASSERT(normalsMatch(mesh, loaded, 1e-4f), "gltf_rt_sphere: normals");
    ASSERT(uvsMatch(mesh, loaded, 1e-4f), "gltf_rt_sphere: UVs");
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
    ASSERT(loaded.hasUVs(), "gltf_rt_sn_lod: UVs");
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
    ASSERT(bromesh::saveGLTF(mesh, path), "gltf_rt_dc_overdraw: save");
    auto scene = bromesh::loadGLTF(path);
    ASSERT(!scene.meshes.empty(), "gltf_rt_dc_overdraw: has meshes");
    auto& loaded = scene.meshes[0];
    ASSERT(loaded.vertexCount() == mesh.vertexCount(), "gltf_rt_dc_overdraw: vertex count");
    ASSERT(loaded.triangleCount() == mesh.triangleCount(), "gltf_rt_dc_overdraw: tri count");
    ASSERT(positionsMatch(mesh, loaded, 1e-4f), "gltf_rt_dc_overdraw: positions");
    ASSERT(loaded.hasUVs(), "gltf_rt_dc_overdraw: UVs preserved");
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
