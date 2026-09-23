#include "test_framework.h"
#include <cmath>
#include <filesystem>
#include <string>
#include <algorithm>
#include <cstring>

#if BROMESH_HAS_GLTF
#include "tiny_gltf.h"

#if BROMESH_HAS_DRACO
#include "bromesh/io/draco.h"
#endif

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

// A CUBICSPLINE sampler's output accessor counts every element, 3 per
// keyframe (in-tangent, value, out-tangent). The saver used to write a third
// of that, and the loader multiplied by 3 again, reading past the accessor.
TEST(gltf_cubicspline_output_count_is_in_elements) {
    auto mesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::Skeleton skeleton;
    bromesh::Bone b0;
    b0.name = "root";
    b0.parent = -1;
    skeleton.bones.push_back(b0);

    bromesh::SkinData skin;
    skin.boneCount = 1;
    skin.boneIndices.assign(mesh.vertexCount() * 4, 0u);
    skin.boneWeights.assign(mesh.vertexCount() * 4, 0.0f);
    for (size_t v = 0; v < mesh.vertexCount(); ++v) skin.boneWeights[v * 4] = 1.0f;
    skin.inverseBindMatrices = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };

    bromesh::Animation anim;
    anim.name = "Cubic";
    anim.duration = 1.0f;
    bromesh::AnimChannel ch;
    ch.boneIndex = 0;
    ch.path = bromesh::AnimChannel::Path::Translation;
    ch.interp = bromesh::AnimChannel::Interp::CubicSpline;
    ch.times = { 0.0f, 1.0f };
    ch.values = {
        0, 0, 0,   1, 2, 3,   0, 0, 0,   // key 0: in, value, out
        0, 0, 0,   4, 5, 6,   0, 0, 0,   // key 1
    };
    anim.channels.push_back(ch);

    std::string path = testDir + "rt_anim_cubic.glb";
    ASSERT(bromesh::saveGLTF(mesh, &skin, &skeleton, {anim}, path), "saveGLTF cubic");

    tinygltf::Model model;
    tinygltf::TinyGLTF reader;
    std::string err, warn;
    ASSERT(reader.LoadBinaryFromFile(&model, &err, &warn, path), "tinygltf reads the saved glb");
    ASSERT(model.animations.size() == 1 && model.animations[0].samplers.size() == 1, "one sampler");
    const auto& s = model.animations[0].samplers[0];
    ASSERT(s.interpolation == "CUBICSPLINE", "sampler is CUBICSPLINE");
    ASSERT(model.accessors[s.input].count == 2, "input count = keyframes");
    ASSERT(model.accessors[s.output].count == 6, "output count = 3 elements per keyframe");

    auto scene = bromesh::loadGLTF(path);
    ASSERT(scene.animations.size() == 1 && scene.animations[0].channels.size() == 1, "loaded channel");
    const auto& lch = scene.animations[0].channels[0];
    ASSERT(lch.values.size() == ch.values.size(), "cubic values not over-read");
    for (size_t i = 0; i < ch.values.size(); ++i)
        ASSERT(std::fabs(lch.values[i] - ch.values[i]) < 1e-6f, "cubic value matches");
    std::remove(path.c_str());
}

namespace {
template <typename T>
void putBytes(std::vector<uint8_t>& buf, size_t at, const T& v) {
    std::memcpy(buf.data() + at, &v, sizeof(T));
}
}

// Every JOINTS_n/WEIGHTS_n set is read and folded into SkinData's 4
// influences: the heaviest 4, duplicate joints summed, renormalized. The
// attributes are stored the awkward ways the spec allows: set 0 interleaved
// with POSITION (bufferView byteStride), set 1 as uint16 joints plus
// normalized uint8 weights.
TEST(gltf_multi_influence_sets_interleaved) {
    const size_t V = 3;
    const size_t stride = 12 + 4 + 16;              // pos f32x3, joints u8x4, weights f32x4
    const size_t interleavedBytes = stride * V;
    const size_t j1Off = interleavedBytes;           // u16x4
    const size_t w1Off = j1Off + 8 * V;              // u8x4 normalized
    std::vector<uint8_t> bytes(w1Off + 4 * V, 0);

    const float pos[V][3] = { {0, 0, 0}, {1, 0, 0}, {0, 0, 1} };
    const uint8_t j0[V][4] = { {0, 1, 2, 3}, {0, 0, 0, 0}, {2, 3, 1, 0} };
    const float w0[V][4] = { {0.3f, 0.1f, 0.05f, 0.05f}, {1, 0, 0, 0}, {0.25f, 0.25f, 0.25f, 0} };
    const uint16_t j1[V][4] = { {4, 5, 6, 7}, {0, 0, 0, 0}, {2, 9, 0, 0} };
    const uint8_t w1[V][4] = { {51, 102, 0, 0}, {0, 0, 0, 0}, {51, 13, 0, 0} };
    for (size_t v = 0; v < V; ++v) {
        for (int c = 0; c < 3; ++c) putBytes(bytes, v * stride + c * 4, pos[v][c]);
        for (int c = 0; c < 4; ++c) bytes[v * stride + 12 + c] = j0[v][c];
        for (int c = 0; c < 4; ++c) putBytes(bytes, v * stride + 16 + c * 4, w0[v][c]);
        for (int c = 0; c < 4; ++c) putBytes(bytes, j1Off + v * 8 + c * 2, j1[v][c]);
        for (int c = 0; c < 4; ++c) bytes[w1Off + v * 4 + c] = w1[v][c];
    }

    tinygltf::Model model;
    model.asset.version = "2.0";
    tinygltf::Buffer buffer;
    buffer.data = bytes;
    model.buffers.push_back(buffer);

    auto addView = [&](size_t off, size_t len, int byteStride) {
        tinygltf::BufferView bv;
        bv.buffer = 0;
        bv.byteOffset = off;
        bv.byteLength = len;
        bv.byteStride = byteStride;
        bv.target = TINYGLTF_TARGET_ARRAY_BUFFER;
        model.bufferViews.push_back(bv);
        return static_cast<int>(model.bufferViews.size()) - 1;
    };
    auto addAcc = [&](int view, size_t off, int compType, int type, bool normalized) {
        tinygltf::Accessor a;
        a.bufferView = view;
        a.byteOffset = off;
        a.componentType = compType;
        a.type = type;
        a.count = V;
        a.normalized = normalized;
        if (type == TINYGLTF_TYPE_VEC3) { a.minValues = {0, 0, 0}; a.maxValues = {1, 0, 1}; }
        model.accessors.push_back(a);
        return static_cast<int>(model.accessors.size()) - 1;
    };
    const int vInter = addView(0, interleavedBytes, static_cast<int>(stride));
    const int vJ1 = addView(j1Off, 8 * V, 0);
    const int vW1 = addView(w1Off, 4 * V, 0);

    tinygltf::Primitive prim;
    prim.mode = TINYGLTF_MODE_TRIANGLES;
    prim.attributes["POSITION"]  = addAcc(vInter, 0,  TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_VEC3, false);
    prim.attributes["JOINTS_0"]  = addAcc(vInter, 12, TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE, TINYGLTF_TYPE_VEC4, false);
    prim.attributes["WEIGHTS_0"] = addAcc(vInter, 16, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_VEC4, false);
    prim.attributes["JOINTS_1"]  = addAcc(vJ1, 0, TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT, TINYGLTF_TYPE_VEC4, false);
    prim.attributes["WEIGHTS_1"] = addAcc(vW1, 0, TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE, TINYGLTF_TYPE_VEC4, true);
    tinygltf::Mesh gm;
    gm.primitives.push_back(prim);
    model.meshes.push_back(gm);
    tinygltf::Node node;
    node.mesh = 0;
    model.nodes.push_back(node);
    tinygltf::Scene sc;
    sc.nodes.push_back(0);
    model.scenes.push_back(sc);
    model.defaultScene = 0;

    std::string path = testDir + "multi_influence.glb";
    tinygltf::TinyGLTF writer;
    ASSERT(writer.WriteGltfSceneToFile(&model, path, true, true, true, true), "wrote glb");

    auto scene = bromesh::loadGLTF(path);
    ASSERT(scene.meshes.size() == 1 && scene.skins.size() == 1, "one mesh + skin");
    const auto& m = scene.meshes[0];
    ASSERT(m.vertexCount() == V, "vertex count");
    // Interleaved POSITION reads each vertex, not a tight run of floats.
    for (size_t v = 0; v < V; ++v)
        for (int c = 0; c < 3; ++c)
            ASSERT(m.positions[v * 3 + c] == pos[v][c], "interleaved position");

    const auto& sk = scene.skins[0];
    ASSERT(sk.boneIndices.size() == V * 4 && sk.boneWeights.size() == V * 4, "4 influences per vertex");
    ASSERT(sk.validate(), "skin validates");

    // v0: heaviest four of {j0 .3, j1 .1, j2 .05, j3 .05, j4 .2, j5 .4}.
    const uint32_t e0j[4] = { 5, 0, 4, 1 };
    const float e0w[4] = { 0.4f, 0.3f, 0.2f, 0.1f };
    for (int k = 0; k < 4; ++k) {
        ASSERT(sk.boneIndices[k] == e0j[k], "v0 joint order");
        ASSERT(std::fabs(sk.boneWeights[k] - e0w[k]) < 1e-3f, "v0 weight");
    }
    // v1: a single joint, still weight 1.
    ASSERT(sk.boneIndices[4] == 0 && std::fabs(sk.boneWeights[4] - 1.0f) < 1e-6f, "v1 joint 0 = 1");
    for (int k = 1; k < 4; ++k) ASSERT(sk.boneWeights[4 + k] == 0.0f, "v1 other weights 0");
    // v2: joint 2 appears in both sets and is summed (.25 + .2).
    ASSERT(sk.boneIndices[8] == 2, "v2 duplicate joint summed to the top");
    float sum = 0.0f;
    for (int k = 0; k < 4; ++k) sum += sk.boneWeights[8 + k];
    ASSERT(std::fabs(sum - 1.0f) < 1e-5f, "v2 renormalized");
    ASSERT(std::fabs(sk.boneWeights[8] - 0.45f / (0.45f + 0.5f + 13.0f / 255.0f)) < 1e-4f, "v2 top weight");

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

TEST(gltf_scene_save_and_load_round_trip) {
    bromesh::GltfScene scene;

    // 2 meshes
    auto boxMesh = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::computeNormals(boxMesh);
    bromesh::projectUVs(boxMesh, bromesh::ProjectionType::Box, 1.0f);
    scene.meshes.push_back(boxMesh);

    auto sphereMesh = bromesh::sphere(0.75f, 16, 12);
    bromesh::computeNormals(sphereMesh);
    bromesh::projectUVs(sphereMesh, bromesh::ProjectionType::Spherical, 1.0f);
    scene.meshes.push_back(sphereMesh);

    // 1 skeleton
    bromesh::Skeleton skel;
    bromesh::Bone b0;
    b0.name = "root";
    b0.parent = -1;
    b0.localT[0] = 0.0f; b0.localT[1] = 0.0f; b0.localT[2] = 0.0f;
    b0.localR[0] = 0.0f; b0.localR[1] = 0.0f; b0.localR[2] = 0.0f; b0.localR[3] = 1.0f;
    b0.localS[0] = 1.0f; b0.localS[1] = 1.0f; b0.localS[2] = 1.0f;
    skel.bones.push_back(b0);
    scene.skeletons.push_back(skel);

    // Skinning data for mesh 0
    bromesh::SkinData skin0;
    skin0.boneCount = 1;
    skin0.boneIndices.assign(boxMesh.vertexCount() * 4, 0u);
    skin0.boneWeights.assign(boxMesh.vertexCount() * 4, 0.0f);
    for (size_t v = 0; v < boxMesh.vertexCount(); ++v) {
        skin0.boneWeights[v * 4 + 0] = 1.0f;
    }
    float idMat[16] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
    skin0.inverseBindMatrices.assign(16, 0.0f);
    std::memcpy(skin0.inverseBindMatrices.data(), idMat, sizeof(idMat));
    scene.skins.push_back(skin0);

    // Empty skin for mesh 1
    scene.skins.push_back({});

    scene.meshSkeleton = { 0, -1 };

    // 1 embedded 2x2 image (RGBA8)
    bromesh::Image img;
    img.name = "test_tex";
    img.width = 2;
    img.height = 2;
    img.mimeType = "image/png";
    img.data = {
        255, 0, 0, 255,     0, 255, 0, 255,
        0, 0, 255, 255,     255, 255, 0, 255
    };
    scene.images.push_back(img);

    // 1 material with baseColorFactor and referencing image 0
    bromesh::Material mat;
    mat.name = "TestPbrMaterial";
    mat.baseColorFactor[0] = 0.8f;
    mat.baseColorFactor[1] = 0.2f;
    mat.baseColorFactor[2] = 0.4f;
    mat.baseColorFactor[3] = 1.0f;
    mat.metallicFactor = 0.5f;
    mat.roughnessFactor = 0.25f;
    mat.emissiveFactor[0] = 0.1f;
    mat.emissiveFactor[1] = 0.2f;
    mat.emissiveFactor[2] = 0.3f;
    mat.baseColorTexture = 0;
    scene.materials.push_back(mat);

    scene.meshMaterial = { 0, -1 };

    // 1 animation targeting skeleton 0
    bromesh::Animation anim;
    anim.name = "SimpleTranslation";
    anim.duration = 1.0f;
    bromesh::AnimChannel ch;
    ch.boneIndex = 0;
    ch.path = bromesh::AnimChannel::Path::Translation;
    ch.interp = bromesh::AnimChannel::Interp::Linear;
    ch.times = { 0.0f, 1.0f };
    ch.values = { 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f };
    anim.channels.push_back(ch);
    scene.animations.push_back(anim);
    scene.animationSkeleton = { 0 };

    std::string path = std::string(testDir) + "rt_scene_multi.glb";
    ASSERT(bromesh::saveGLTF(scene, path), "saveGLTF(GltfScene) returns true");

    auto loaded = bromesh::loadGLTF(path);
    ASSERT(loaded.meshes.size() == 2, "2 meshes loaded");
    ASSERT(positionsMatch(boxMesh, loaded.meshes[0], 1e-4f), "mesh 0 positions match");
    ASSERT(positionsMatch(sphereMesh, loaded.meshes[1], 1e-4f), "mesh 1 positions match");
    ASSERT(loaded.meshes[0].hasNormals() && loaded.meshes[0].hasUVs(), "mesh 0 has normals and UVs");
    ASSERT(loaded.meshes[1].hasNormals() && loaded.meshes[1].hasUVs(), "mesh 1 has normals and UVs");

    // Skeletons and skins
    ASSERT(loaded.skeletons.size() == 1, "1 skeleton loaded");
    ASSERT(loaded.meshSkeleton.size() == 2, "meshSkeleton size 2");
    ASSERT(loaded.meshSkeleton[0] == 0, "mesh 0 has skeleton 0");
    ASSERT(loaded.meshSkeleton[1] == -1, "mesh 1 has skeleton -1");
    ASSERT(loaded.skins.size() == 2, "skins size 2");
    ASSERT(loaded.skins[0].boneIndices.size() == boxMesh.vertexCount() * 4, "mesh 0 skin indices size");

    // Materials
    ASSERT(loaded.materials.size() == 1, "1 material loaded");
    ASSERT(loaded.materials[0].name == "TestPbrMaterial", "material name matches");
    ASSERT(std::fabs(loaded.materials[0].baseColorFactor[0] - 0.8f) < 1e-3f, "baseColorFactor R");
    ASSERT(std::fabs(loaded.materials[0].baseColorFactor[1] - 0.2f) < 1e-3f, "baseColorFactor G");
    ASSERT(std::fabs(loaded.materials[0].baseColorFactor[2] - 0.4f) < 1e-3f, "baseColorFactor B");
    ASSERT(std::fabs(loaded.materials[0].metallicFactor - 0.5f) < 1e-3f, "metallicFactor");
    ASSERT(std::fabs(loaded.materials[0].roughnessFactor - 0.25f) < 1e-3f, "roughnessFactor");
    ASSERT(loaded.materials[0].baseColorTexture == 0, "baseColorTexture index resolved to image 0");

    // Mesh material assignment
    ASSERT(loaded.meshMaterial.size() == 2, "meshMaterial size 2");
    ASSERT(loaded.meshMaterial[0] == 0, "mesh 0 assigned material 0");
    ASSERT(loaded.meshMaterial[1] == -1, "mesh 1 assigned material -1");

    // Images
    ASSERT(loaded.images.size() == 1, "1 image loaded");
    ASSERT(loaded.images[0].width == 2 && loaded.images[0].height == 2, "image dimensions 2x2");
    ASSERT(loaded.images[0].data.size() == 16, "image data size 16 (2x2x4)");
    ASSERT(loaded.images[0].data == img.data, "embedded image pixels match exactly");

    // Animations
    ASSERT(loaded.animations.size() == 1, "1 animation loaded");
    ASSERT(loaded.animations[0].name == "SimpleTranslation", "animation name matches");
    ASSERT(loaded.animations[0].channels.size() == 1, "1 anim channel");
    ASSERT(loaded.animations[0].channels[0].times.size() == 2, "channel times count");
    ASSERT(loaded.animations[0].channels[0].values.size() == 6, "channel values count");

    std::remove(path.c_str());
}

#if BROMESH_HAS_DRACO && BROMESH_HAS_GLTF
TEST(gltf_draco_extension_loading) {
    auto box = bromesh::box(1.0f, 0.8f, 0.6f);
    bromesh::computeNormals(box);
    bromesh::projectUVs(box, bromesh::ProjectionType::Box, 1.0f);
    ASSERT(box.validate(), "box validates");

    std::string err;
    std::vector<uint8_t> dracoBytes = bromesh::encodeDraco(box, bromesh::DracoEncodeOptions(), &err);
    ASSERT(!dracoBytes.empty(), "draco encode succeeded");

    // Build tinygltf::Model using KHR_draco_mesh_compression
    tinygltf::Model model;
    model.asset.version = "2.0";
    model.asset.generator = "bromesh_draco_test";

    // Buffer 0 holds dracoBytes
    tinygltf::Buffer buffer;
    buffer.data = dracoBytes;
    model.buffers.push_back(buffer);

    // BufferView 0 for the draco compressed payload
    tinygltf::BufferView bv;
    bv.buffer = 0;
    bv.byteOffset = 0;
    bv.byteLength = dracoBytes.size();
    model.bufferViews.push_back(bv);

    // In glTF, standard accessors are provided as fallbacks for loaders without Draco
    // For Draco-enabled loaders, accessors specify count and component type
    tinygltf::Accessor posAcc;
    posAcc.bufferView = -1;
    posAcc.count = box.vertexCount();
    posAcc.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
    posAcc.type = TINYGLTF_TYPE_VEC3;
    model.accessors.push_back(posAcc);

    tinygltf::Accessor nrmAcc;
    nrmAcc.bufferView = -1;
    nrmAcc.count = box.vertexCount();
    nrmAcc.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
    nrmAcc.type = TINYGLTF_TYPE_VEC3;
    model.accessors.push_back(nrmAcc);

    tinygltf::Accessor uvAcc;
    uvAcc.bufferView = -1;
    uvAcc.count = box.vertexCount();
    uvAcc.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
    uvAcc.type = TINYGLTF_TYPE_VEC2;
    model.accessors.push_back(uvAcc);

    tinygltf::Accessor idxAcc;
    idxAcc.bufferView = -1;
    idxAcc.count = box.indices.size();
    idxAcc.componentType = TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT;
    idxAcc.type = TINYGLTF_TYPE_SCALAR;
    model.accessors.push_back(idxAcc);

    tinygltf::Mesh mesh;
    tinygltf::Primitive prim;
    prim.mode = TINYGLTF_MODE_TRIANGLES;
    prim.indices = 3;
    prim.attributes["POSITION"] = 0;
    prim.attributes["NORMAL"] = 1;
    prim.attributes["TEXCOORD_0"] = 2;

    // Draco extension dictionary: bufferView = 0, attributes dictionary
    tinygltf::Value::Object dracoExtObj;
    dracoExtObj["bufferView"] = tinygltf::Value(0);
    tinygltf::Value::Object attrsObj;
    attrsObj["POSITION"] = tinygltf::Value(0);
    attrsObj["NORMAL"] = tinygltf::Value(1);
    attrsObj["TEXCOORD_0"] = tinygltf::Value(2);
    dracoExtObj["attributes"] = tinygltf::Value(attrsObj);

    prim.extensions["KHR_draco_mesh_compression"] = tinygltf::Value(dracoExtObj);
    mesh.primitives.push_back(prim);
    model.meshes.push_back(mesh);

    tinygltf::Node node;
    node.mesh = 0;
    model.nodes.push_back(node);

    tinygltf::Scene sceneNode;
    sceneNode.nodes.push_back(0);
    model.scenes.push_back(sceneNode);
    model.defaultScene = 0;

    model.extensionsUsed.push_back("KHR_draco_mesh_compression");
    model.extensionsRequired.push_back("KHR_draco_mesh_compression");

    std::string path = std::string(testDir) + "draco_ext_test.glb";
    tinygltf::TinyGLTF writer;
    ASSERT(writer.WriteGltfSceneToFile(&model, path, true, true, true, true), "wrote draco glb");

    // Load with loadGLTF
    auto loaded = bromesh::loadGLTF(path);
    ASSERT(loaded.meshes.size() == 1, "loaded 1 mesh from draco compressed glb");
    const auto& loadedMesh = loaded.meshes[0];
    ASSERT(loadedMesh.triangleCount() == box.triangleCount(), "triangle count matches");
    ASSERT(loadedMesh.hasNormals(), "loaded draco mesh has normals");
    ASSERT(loadedMesh.hasUVs(), "loaded draco mesh has UVs");
    ASSERT(loadedMesh.validate(), "loaded draco mesh validates");

    // Extents check (Draco quantization nudges positions slightly)
    float lo[3] = {1e9f, 1e9f, 1e9f}, hi[3] = {-1e9f, -1e9f, -1e9f};
    float dlo[3] = {1e9f, 1e9f, 1e9f}, dhi[3] = {-1e9f, -1e9f, -1e9f};
    for (size_t i = 0; i < box.positions.size(); i += 3) {
        for (int c = 0; c < 3; ++c) {
            lo[c] = std::min(lo[c], box.positions[i + c]);
            hi[c] = std::max(hi[c], box.positions[i + c]);
        }
    }
    for (size_t i = 0; i < loadedMesh.positions.size(); i += 3) {
        for (int c = 0; c < 3; ++c) {
            dlo[c] = std::min(dlo[c], loadedMesh.positions[i + c]);
            dhi[c] = std::max(dhi[c], loadedMesh.positions[i + c]);
        }
    }
    for (int c = 0; c < 3; ++c) {
        ASSERT(std::fabs(lo[c] - dlo[c]) < 1e-3f, "min bounds survive draco compression");
        ASSERT(std::fabs(hi[c] - dhi[c]) < 1e-3f, "max bounds survive draco compression");
    }

    std::remove(path.c_str());
}
#endif

#endif // BROMESH_HAS_GLTF

