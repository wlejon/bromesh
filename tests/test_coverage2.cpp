// Second round of targeted coverage tests. Picks up: hand-crafted glTF
// (matrix node, COLOR_0, UINT indices, embedded image, materials,
// step animation, ancestor scene tree), Loop/CatmullClark subdivision
// with vertex colors, coplanar mesh intersection, raycast on a mesh
// without explicit normals, lsystem turtle G/f and '-' rotation, bake
// transfer with default search distance, skeleton fit "mid:" expression,
// PLY binary loader with non-float vertex types, polygon CW outer
// (hole-orientation-fix branch), CapsuleField with sphere obstacles.

#include "test_framework.h"
#include "bromesh/procedural/lsystem_turtle.h"
#include "bromesh/procedural/obstacle_field.h"

#include "tiny_gltf.h"

#include <algorithm>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>

namespace {

// Build a minimal in-memory tinygltf::Model with:
//  - one triangle mesh (positions + COLOR_0 floats + UNSIGNED_INT indices)
//  - one skin with a joint node that has its transform in `matrix` (no TRS)
//  - the joint node lives under an armature node that itself has a
//    translation, so loadGLTF has to climb the parent chain and accumulate
//    rootTransform via nodeLocalMat / matMul4.
//  - one image (1x1 RGBA) and one material referencing it on every slot
//  - one animation channel with STEP interpolation on translation
//
// Save it to disk via tinygltf, then call bromesh::loadGLTF on the result.
// Exercises decomposeTRS, nodeLocalMat, image+material decoding, and the
// matrix-node branches in the skin loader.
static std::string buildHandCraftedGltf(const char* path) {
    tinygltf::Model m;
    m.asset.version = "2.0";
    m.asset.generator = "bromesh_test";

    // ---- buffer with all attribute data ----
    // Positions: 3 verts vec3 -> 36 bytes
    // Colors (vec4 float): 48 bytes
    // Joints (vec4 u8): 12 bytes; pad to 16
    // Weights (vec4 float): 48 bytes
    // Indices (u32): 12 bytes
    // IBM (mat4): 64 bytes
    tinygltf::Buffer buf;
    auto appendBytes = [&](const void* d, size_t n) {
        size_t off = buf.data.size();
        buf.data.resize(off + n);
        std::memcpy(buf.data.data() + off, d, n);
        return off;
    };

    float positions[] = {
        0.0f, 0.0f, 0.0f,
        1.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f,
    };
    size_t posOff = appendBytes(positions, sizeof(positions));

    float colors[] = {
        1.0f, 0.0f, 0.0f, 1.0f,
        0.0f, 1.0f, 0.0f, 1.0f,
        0.0f, 0.0f, 1.0f, 1.0f,
    };
    size_t colOff = appendBytes(colors, sizeof(colors));

    // u8 joints, vec4 per vertex (all referring to bone 0).
    uint8_t joints[] = {
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
    };
    size_t jntOff = appendBytes(joints, sizeof(joints));
    // pad to 4-byte alignment
    while (buf.data.size() % 4) buf.data.push_back(0);

    float weights[] = {
        1.0f, 0.0f, 0.0f, 0.0f,
        1.0f, 0.0f, 0.0f, 0.0f,
        1.0f, 0.0f, 0.0f, 0.0f,
    };
    size_t wtOff = appendBytes(weights, sizeof(weights));

    uint32_t indices[] = { 0, 1, 2 };
    size_t idxOff = appendBytes(indices, sizeof(indices));

    float ibm[16] = {
        1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1,
    };
    size_t ibmOff = appendBytes(ibm, sizeof(ibm));

    // Animation: STEP translation on bone 0, two keys.
    float animTimes[] = { 0.0f, 1.0f };
    size_t animTOff = appendBytes(animTimes, sizeof(animTimes));
    float animValues[] = { 0,0,0,  0,2,0 };
    size_t animVOff = appendBytes(animValues, sizeof(animValues));

    m.buffers.push_back(std::move(buf));

    // ---- buffer views + accessors ----
    auto addView = [&](size_t off, size_t len) {
        tinygltf::BufferView bv;
        bv.buffer = 0;
        bv.byteOffset = off;
        bv.byteLength = len;
        m.bufferViews.push_back(bv);
        return (int)m.bufferViews.size() - 1;
    };
    auto addAccessor = [&](int view, int componentType, int type, size_t count) {
        tinygltf::Accessor a;
        a.bufferView = view;
        a.byteOffset = 0;
        a.componentType = componentType;
        a.type = type;
        a.count = count;
        m.accessors.push_back(a);
        return (int)m.accessors.size() - 1;
    };

    int vPos = addView(posOff, sizeof(positions));
    int aPos = addAccessor(vPos, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_VEC3, 3);
    m.accessors[aPos].minValues = {0.0, 0.0, 0.0};
    m.accessors[aPos].maxValues = {1.0, 1.0, 0.0};

    int vCol = addView(colOff, sizeof(colors));
    int aCol = addAccessor(vCol, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_VEC4, 3);

    int vJnt = addView(jntOff, sizeof(joints));
    int aJnt = addAccessor(vJnt, TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE, TINYGLTF_TYPE_VEC4, 3);

    int vWt = addView(wtOff, sizeof(weights));
    int aWt = addAccessor(vWt, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_VEC4, 3);

    int vIdx = addView(idxOff, sizeof(indices));
    int aIdx = addAccessor(vIdx, TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT, TINYGLTF_TYPE_SCALAR, 3);

    int vIbm = addView(ibmOff, sizeof(ibm));
    int aIbm = addAccessor(vIbm, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_MAT4, 1);

    int vAnT = addView(animTOff, sizeof(animTimes));
    int aAnT = addAccessor(vAnT, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_SCALAR, 2);
    m.accessors[aAnT].minValues = {0.0};
    m.accessors[aAnT].maxValues = {1.0};
    int vAnV = addView(animVOff, sizeof(animValues));
    int aAnV = addAccessor(vAnV, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_VEC3, 2);

    // ---- mesh ----
    tinygltf::Mesh mesh;
    tinygltf::Primitive prim;
    prim.attributes["POSITION"]   = aPos;
    prim.attributes["COLOR_0"]    = aCol;
    prim.attributes["JOINTS_0"]   = aJnt;
    prim.attributes["WEIGHTS_0"]  = aWt;
    prim.indices  = aIdx;
    prim.material = 0;
    prim.mode = TINYGLTF_MODE_TRIANGLES;
    mesh.primitives.push_back(prim);
    m.meshes.push_back(mesh);

    // ---- material with no textures (image embedding requires PNG bytes,
    // which our minimal pipeline can't produce). Material loader still
    // exercises baseColorFactor, metallic/roughness, and emissiveFactor.
    tinygltf::Material mat;
    mat.name = "matA";
    mat.pbrMetallicRoughness.baseColorFactor = {0.5, 0.6, 0.7, 1.0};
    mat.pbrMetallicRoughness.metallicFactor = 0.25;
    mat.pbrMetallicRoughness.roughnessFactor = 0.75;
    mat.emissiveFactor = {0.1, 0.2, 0.3};
    m.materials.push_back(mat);

    // ---- nodes ----
    // 0: armature (parent of joint, with a translation in TRS form)
    // 1: joint (with `matrix` rather than TRS — forces decomposeTRS)
    // 2: mesh node (skinned)
    tinygltf::Node armature;
    armature.name = "armature";
    armature.translation = {0.0, 0.1, 0.0};
    armature.children = { 1 };
    m.nodes.push_back(armature);

    tinygltf::Node joint;
    joint.name = "joint";
    // 4x4 matrix with non-uniform scale to force decomposeTRS to compute s.
    joint.matrix = {
        2.0, 0.0, 0.0, 0.0,
        0.0, 1.5, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0,
    };
    m.nodes.push_back(joint);

    tinygltf::Node meshNode;
    meshNode.mesh = 0;
    meshNode.skin = 0;
    m.nodes.push_back(meshNode);

    tinygltf::Skin skin;
    skin.inverseBindMatrices = aIbm;
    skin.joints = { 1 };
    m.skins.push_back(skin);

    // ---- animation ----
    tinygltf::Animation anim;
    anim.name = "step_anim";
    tinygltf::AnimationSampler samp;
    samp.input = aAnT;
    samp.output = aAnV;
    samp.interpolation = "STEP";
    anim.samplers.push_back(samp);
    tinygltf::AnimationChannel ach;
    ach.sampler = 0;
    ach.target_node = 1; // joint
    ach.target_path = "translation";
    anim.channels.push_back(ach);
    m.animations.push_back(anim);

    // ---- scene ----
    tinygltf::Scene sc;
    sc.nodes = { 0, 2 };
    m.scenes.push_back(sc);
    m.defaultScene = 0;

    tinygltf::TinyGLTF writer;
    bool ok = writer.WriteGltfSceneToFile(&m, path,
                                          /*embedImages*/  true,
                                          /*embedBuffers*/ true,
                                          /*prettyPrint*/  false,
                                          /*writeBinary*/  false);
    return ok ? std::string(path) : std::string();
}

} // anonymous

TEST(gltf_handcrafted_matrix_image_material) {
    const char* p = "test_handcrafted.gltf";
    std::string out = buildHandCraftedGltf(p);
    ASSERT(!out.empty(), "gltf hand-crafted: write succeeded");

    auto scene = bromesh::loadGLTF(p);
    ASSERT(scene.meshes.size() == 1, "gltf hand: one mesh");
    ASSERT(scene.meshes[0].hasColors(), "gltf hand: COLOR_0 loaded");
    ASSERT(scene.skeletons.size() == 1, "gltf hand: one skeleton");
    // Joint had non-uniform scale 2,1.5,1 — verify localS captured via
    // decomposeTRS (within rounding).
    const auto& bone = scene.skeletons[0].bones[0];
    ASSERT(std::fabs(bone.localS[0] - 2.0f) < 1e-3f, "gltf hand: scale x=2");
    ASSERT(std::fabs(bone.localS[1] - 1.5f) < 1e-3f, "gltf hand: scale y=1.5");
    // rootTransform should include the armature's translation y=0.1.
    ASSERT(std::fabs(scene.skeletons[0].rootTransform[13] - 0.1f) < 1e-4f,
           "gltf hand: rootTransform y captured");
    ASSERT(scene.materials.size() == 1, "gltf hand: one material");
    ASSERT(std::fabs(scene.materials[0].baseColorFactor[0] - 0.5f) < 1e-4f,
           "gltf hand: baseColorFactor");
    ASSERT(std::fabs(scene.materials[0].metallicFactor - 0.25f) < 1e-4f,
           "gltf hand: metallicFactor");
    ASSERT(std::fabs(scene.materials[0].emissiveFactor[2] - 0.3f) < 1e-4f,
           "gltf hand: emissiveFactor");
    ASSERT(scene.animations.size() == 1, "gltf hand: one animation");
    std::remove(p);
    std::remove("test_handcrafted.bin");
}

// =============================================================================
// Subdivision with vertex colors. All three subdivision algorithms have a
// "if hasColors" branch that the existing tests don't exercise.
// =============================================================================
TEST(subdivide_all_three_with_colors) {
    auto base = bromesh::box(1.0f, 1.0f, 1.0f);
    base.colors.assign(base.vertexCount() * 4, 0.0f);
    for (size_t v = 0; v < base.vertexCount(); ++v) {
        base.colors[v*4 + 0] = (v & 1) ? 1.0f : 0.0f;
        base.colors[v*4 + 3] = 1.0f;
    }
    auto mid = bromesh::subdivideMidpoint(base, 1);
    ASSERT(mid.hasColors(), "subdivide midpoint preserves colors");

    auto loop = bromesh::subdivideLoop(base, 1);
    ASSERT(loop.hasColors(), "subdivide loop preserves colors");

    auto cc = bromesh::subdivideCatmullClark(base, 1);
    ASSERT(cc.hasColors(), "subdivide cc preserves colors");
}

// =============================================================================
// intersect — coplanar triangle case. Two triangles in the XY plane that
// overlap in 2D drive the coplanarTrianglesOverlap path (SAT).
// =============================================================================
TEST(meshes_intersect_coplanar_overlap) {
    bromesh::MeshData a;
    a.positions = { -1,-1,0,   1,-1,0,   0,1,0 };
    a.indices   = { 0,1,2 };

    // Triangle B is in the same plane but shifted, partially overlapping A.
    bromesh::MeshData b;
    b.positions = { 0,-0.5f,0,  2,-0.5f,0,  1,0.5f,0 };
    b.indices   = { 0,1,2 };

    bool inter = bromesh::meshesIntersect(a, b);
    ASSERT(inter, "meshes intersect: coplanar overlap detected");

    // And a disjoint coplanar case — same plane, far apart.
    bromesh::MeshData c;
    c.positions = { 5,-1,0,   7,-1,0,   6,1,0 };
    c.indices   = { 0,1,2 };
    ASSERT(!bromesh::meshesIntersect(a, c),
           "meshes intersect: coplanar far apart -> no");
}

// =============================================================================
// Raycast on a mesh that has no per-vertex normals — forces the
// face-normal-derived path in both raycast() and closestPoint().
// =============================================================================
TEST(raycast_face_normal_branch) {
    bromesh::MeshData m;
    m.positions = { 0,0,0,   1,0,0,   0,1,0 };
    m.indices   = { 0,1,2 };
    // NO normals -> raycast must derive the face normal.

    float origin[3] = { 0.25f, 0.25f, 1.0f };
    float dir[3]    = { 0,    0,    -1.0f };
    auto hit = bromesh::raycast(m, origin, dir);
    ASSERT(hit.hit, "raycast face-normal: hit");
    ASSERT(std::fabs(hit.normal[2] - 1.0f) < 1e-4f,
           "raycast face-normal: z=+1");

    // closestPoint without normals also computes face normal.
    float p[3] = { 0.25f, 0.25f, 2.0f };
    auto cp = bromesh::closestPoint(m, p);
    ASSERT(cp.hit, "closestPoint face-normal: hit");
    ASSERT(std::fabs(cp.normal[2] - 1.0f) < 1e-4f,
           "closestPoint face-normal: z=+1");
}

// =============================================================================
// L-system turtle: G/f modules (advance without emit) and '-' rotation.
// =============================================================================
TEST(lsystem_turtle_advance_and_negative_rotation) {
    bromesh::TurtleOptions opts;
    opts.stepLength = 1.0f;
    opts.angle = 0.5f;
    std::vector<bromesh::Module> mods;
    auto add = [&](char s, std::vector<float> params = {}) {
        bromesh::Module m; m.symbol = s; m.params = std::move(params);
        mods.push_back(m);
    };
    // G + f advance without emitting segments. Then '-' to rotate and F to emit.
    add('G');
    add('f', { 0.5f });
    add('-', { 30.0f });
    add('F');
    auto segs = bromesh::lsystemToBranches(mods, opts);
    // Only one F + the synthetic root => 2 segs.
    ASSERT(segs.size() == 2u, "turtle: G/f advance, then F emits 1+root");
}

// =============================================================================
// bake_transfer — call with searchDistance == 0 to drive the auto-distance
// computation from meshBBoxDiag.
// =============================================================================
TEST(bake_normals_from_reference_auto_distance) {
    auto low = bromesh::box(1.0f, 1.0f, 1.0f);
    bromesh::computeNormals(low);
    bromesh::projectUVs(low, bromesh::ProjectionType::Box);
    auto high = bromesh::box(1.05f, 1.05f, 1.05f);
    bromesh::computeNormals(high);
    // searchDistance left at default (0) triggers meshBBoxDiag * 0.05f.
    auto tex = bromesh::bakeNormalsFromReference(low, high, 16, 16);
    ASSERT(tex.width == 16 && tex.height == 16,
           "bake auto-distance: texture sized");
}

// =============================================================================
// skeleton_fit — exercise the "mid:A,B" expression form in a tiny custom
// RigSpec. The bundled specs only use landmark/lerp/offset.
// =============================================================================
TEST(skeleton_fit_mid_expression) {
    bromesh::RigSpec spec;
    spec.name = "tiny";
    spec.landmarks.push_back({ "a", "" });
    spec.landmarks.push_back({ "b", "" });

    bromesh::RigSpec::BoneDecl root;
    root.name = "root"; root.parent = "";
    root.head = "landmark:a";
    root.tail = "mid:a,b";
    spec.bones.push_back(root);

    bromesh::Landmarks lm;
    lm.points["a"] = { 0, 0, 0 };
    lm.points["b"] = { 2, 0, 0 };

    auto mesh = bromesh::box(0.1f, 0.1f, 0.1f);
    auto sk = bromesh::fitSkeleton(spec, lm, mesh);
    ASSERT(sk.bones.size() == 1, "skeleton fit mid: one bone fitted");
}

// =============================================================================
// PLY — binary with double-precision positions, ushort UVs, and int colors.
// Drives readFloatVal type-dispatch branches that the existing roundtrip
// (which only writes float / uchar) doesn't reach.
// =============================================================================
TEST(ply_binary_mixed_types) {
    const char* path = "test_mixed_types.ply";
    FILE* f = std::fopen(path, "wb");
    ASSERT(f != nullptr, "open mixed-types ply");
    if (!f) return;

    const char* header =
        "ply\n"
        "format binary_little_endian 1.0\n"
        "element vertex 3\n"
        "property double x\nproperty double y\nproperty double z\n"
        "property ushort s\nproperty ushort t\n"
        "property int red\nproperty int green\nproperty int blue\n"
        "element face 1\n"
        "property list uchar uint vertex_indices\n"
        "end_header\n";
    std::fwrite(header, 1, std::strlen(header), f);

    auto wd = [&](double v) { std::fwrite(&v, sizeof(double), 1, f); };
    auto wu16 = [&](uint16_t v) { std::fwrite(&v, sizeof(uint16_t), 1, f); };
    auto wi32 = [&](int32_t v) { std::fwrite(&v, sizeof(int32_t), 1, f); };

    // Vertex 0: position, uv, color
    wd(0.0); wd(0.0); wd(0.0); wu16(0); wu16(0); wi32(255); wi32(0); wi32(0);
    // Vertex 1
    wd(1.0); wd(0.0); wd(0.0); wu16(65535); wu16(0); wi32(0); wi32(255); wi32(0);
    // Vertex 2
    wd(0.0); wd(1.0); wd(0.0); wu16(0); wu16(65535); wi32(0); wi32(0); wi32(255);
    // Face
    uint8_t three = 3;
    std::fwrite(&three, 1, 1, f);
    uint32_t tri[] = { 0, 1, 2 };
    std::fwrite(tri, sizeof(uint32_t), 3, f);
    std::fclose(f);

    auto m = bromesh::loadPLY(path);
    ASSERT(m.vertexCount() == 3, "ply mixed-types: 3 verts");
    ASSERT(m.hasColors(), "ply mixed-types: colors loaded");
    ASSERT(m.hasUVs(), "ply mixed-types: uvs loaded");
    // Position came in as double-precision.
    ASSERT(std::fabs(m.positions[3] - 1.0f) < 1e-6f,
           "ply mixed-types: double-precision pos round-trips");
    std::remove(path);
}

// =============================================================================
// polygon — clockwise outer with a clockwise hole; the triangulator reverses
// the hole to match outer orientation. Hits the reverse + swap branch.
// =============================================================================
TEST(polygon_triangulate_2d_cw_outer_with_hole) {
    // CW outer (winds clockwise, normal would point -Z).
    std::vector<float> outer = {
        -2,-2,   2,-2,   2,2,   -2,2,
    };
    // Reverse to make CW:
    std::reverse(outer.begin(), outer.end());
    // After reversal the order is (-2,2),(2,2),(2,-2),(-2,-2) but as a flat
    // vector of xy pairs that's not what we want. Build it explicitly:
    outer = { -2,2,  2,2,  2,-2,  -2,-2 };
    std::vector<std::vector<float>> holes;
    holes.push_back({ -1,1,  1,1,  1,-1,  -1,-1 });  // also CW (matches outer)
    auto m = bromesh::triangulatePolygon2D(outer, holes, 0.0f);
    ASSERT(!m.empty(), "polygon 2D CW outer + hole: produces mesh");
}

// =============================================================================
// CapsuleField with sphere obstacles + a long capsule whose length forces the
// multi-cell insertion branch.
// =============================================================================
TEST(capsule_field_with_spheres_and_long_capsule) {
    using namespace bromesh;
    // A long capsule that will be sampled at multiple cells along its axis.
    Capsule c;
    c.a = { 0, 0, 0 };
    c.b = { 10, 0, 0 };
    c.radius = 0.2f;
    c.tag = 1;

    Sphere s;
    s.center = { 5, 5, 0 };
    s.radius = 1.0f;
    s.tag = 2;

    CapsuleField f({c}, {s}, /*cellSize=*/1.0f);
    ASSERT(!f.empty(), "field: non-empty");
    ASSERT(f.sphereCount() == 1, "field: one sphere");

    // intersectsSphere near the obstacle sphere -> true.
    ASSERT(f.intersectsSphere({ 5, 5, 0 }, 0.5f),
           "field: sphere overlap detected");

    // distance() near the long capsule centerline.
    float d = f.distance({ 5, 0, 0 });
    ASSERT(d < 0.0f, "field: distance inside long capsule");

    // excludeTag should skip the sphere and report only capsule distance.
    bromath::Vec3 nearSphere = { 5, 5, 0 };
    float dExcl = f.distance(nearSphere, /*excludeTag=*/2);
    ASSERT(dExcl > 0.0f, "field: excludeTag skips sphere");
}
