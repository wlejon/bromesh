#include "test_framework.h"

// Rig / animation / transfer / bake tests (new additions)
// =============================================================================

// Build a trivial 2-bone skeleton: root at origin, child 1 unit up.
static bromesh::Skeleton makeTwoBoneSkeleton() {
    bromesh::Skeleton sk;
    bromesh::Bone root;
    root.name = "root"; root.parent = -1;
    root.localT[0]=0; root.localT[1]=0; root.localT[2]=0;
    // Identity inverse bind
    sk.bones.push_back(root);

    bromesh::Bone child;
    child.name = "child"; child.parent = 0;
    child.localT[0]=0; child.localT[1]=1; child.localT[2]=0;
    // Inverse bind translates so child is at origin in its own space
    child.inverseBind[12] = 0; child.inverseBind[13] = -1; child.inverseBind[14] = 0;
    sk.bones.push_back(child);
    return sk;
}

TEST(skeleton_find_bone) {
    auto sk = makeTwoBoneSkeleton();
    ASSERT(sk.findBone("root") == 0, "findBone root");
    ASSERT(sk.findBone("child") == 1, "findBone child");
    ASSERT(sk.findBone("nope") == -1, "findBone missing");
}

TEST(gltf_skin_roundtrip) {
    auto sk = makeTwoBoneSkeleton();
    auto mesh = bromesh::box(0.5f, 0.5f, 0.5f);

    bromesh::SkinData skin;
    skin.boneCount = 2;
    skin.boneIndices.assign(mesh.vertexCount() * 4, 0u);
    skin.boneWeights.assign(mesh.vertexCount() * 4, 0.0f);
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        // Upper verts bind to child, lower to root
        if (mesh.positions[v*3+1] > 0) skin.boneIndices[v*4+0] = 1;
        skin.boneWeights[v*4+0] = 1.0f;
    }

    bromesh::Animation anim;
    anim.name = "wiggle";
    bromesh::AnimChannel ch;
    ch.boneIndex = 1;
    ch.path = bromesh::AnimChannel::Path::Rotation;
    ch.interp = bromesh::AnimChannel::Interp::Linear;
    ch.times = { 0.0f, 1.0f };
    ch.values = { 0,0,0,1,  0,0.7071f,0,0.7071f };
    anim.channels.push_back(ch);
    anim.duration = 1.0f;
    std::vector<bromesh::Animation> anims = { anim };

    const char* tmpPath = "test_skinned.glb";
    bool saved = bromesh::saveGLTF(mesh, &skin, &sk, anims, tmpPath);
    ASSERT(saved, "saveGLTF with skin/skel/anim should succeed");

    auto loaded = bromesh::loadGLTF(tmpPath);
    ASSERT(loaded.meshes.size() == 1, "one mesh loaded");
    ASSERT(loaded.skeletons.size() == 1, "one skeleton loaded");
    ASSERT(loaded.skeletons[0].bones.size() == 2, "two bones loaded");
    ASSERT(loaded.animations.size() == 1, "one animation loaded");
    ASSERT(loaded.animations[0].channels.size() == 1, "one channel loaded");
    std::remove(tmpPath);
}

TEST(pose_bind_identity) {
    auto sk = makeTwoBoneSkeleton();
    auto pose = bromesh::bindPose(sk);
    std::vector<float> skin;
    bromesh::computeSkinningMatrices(sk, pose, skin);
    ASSERT(skin.size() == 32, "2 bones * 16 floats");
    // With matching bind/inverse-bind, skinning matrices should be identity
    ASSERT(std::fabs(skin[0] - 1.0f) < 1e-4f, "skin[0][0] == 1 at bind");
    ASSERT(std::fabs(skin[5] - 1.0f) < 1e-4f, "skin[1][1] == 1 at bind");
    ASSERT(std::fabs(skin[10] - 1.0f) < 1e-4f, "skin[2][2] == 1 at bind");
}

TEST(pose_evaluate_linear) {
    auto sk = makeTwoBoneSkeleton();
    bromesh::Animation anim;
    anim.name = "t";
    bromesh::AnimChannel ch;
    ch.boneIndex = 1;
    ch.path = bromesh::AnimChannel::Path::Translation;
    ch.interp = bromesh::AnimChannel::Interp::Linear;
    ch.times = { 0.0f, 1.0f };
    ch.values = { 0,1,0,  0,3,0 };
    anim.channels.push_back(ch);
    anim.duration = 1.0f;

    auto pose = bromesh::evaluateAnimation(sk, anim, 0.5f, false);
    // Expect child local translation y == 2 at t=0.5
    ASSERT(std::fabs(pose.data[1*10 + 1] - 2.0f) < 1e-3f, "linear lerp midpoint");
}

TEST(weight_transfer_identity) {
    auto sk = makeTwoBoneSkeleton();
    // Source: a simple plane (2 tris), skinned: left->root, right->child
    bromesh::MeshData src;
    src.positions = { -1,0,0,  1,0,0,  -1,2,0,  1,2,0 };
    src.normals = { 0,0,1, 0,0,1, 0,0,1, 0,0,1 };
    src.indices = { 0, 1, 2,  2, 1, 3 };

    bromesh::SkinData srcSkin;
    srcSkin.boneCount = 2;
    srcSkin.boneIndices = { 0,0,0,0,  0,0,0,0,  1,0,0,0,  1,0,0,0 };
    srcSkin.boneWeights = { 1,0,0,0,  1,0,0,0,  1,0,0,0,  1,0,0,0 };

    // Target: two points, one on each side
    bromesh::MeshData tgt;
    tgt.positions = { -0.5f, 0.01f, 0,  -0.5f, 1.9f, 0 };
    tgt.indices = { 0, 1, 0 }; // degenerate but present

    auto out = bromesh::transferSkinWeights(tgt, src, srcSkin);
    ASSERT(out.boneIndices.size() == 8, "4 weights * 2 verts");
    // First vertex (low Y) should bind to bone 0
    ASSERT(out.boneIndices[0*4+0] == 0, "low vertex -> root");
    // Second vertex (high Y) should bind to bone 1
    ASSERT(out.boneIndices[1*4+0] == 1, "high vertex -> child");
    // Weights sum to 1
    float s0 = out.boneWeights[0] + out.boneWeights[1] + out.boneWeights[2] + out.boneWeights[3];
    ASSERT(std::fabs(s0 - 1.0f) < 1e-3f, "weights sum to 1 v0");
    float s1 = out.boneWeights[4] + out.boneWeights[5] + out.boneWeights[6] + out.boneWeights[7];
    ASSERT(std::fabs(s1 - 1.0f) < 1e-3f, "weights sum to 1 v1");
}

TEST(shrinkwrap_to_plane) {
    bromesh::MeshData plane;
    plane.positions = { -2, 0, -2,  2, 0, -2,  -2, 0, 2,  2, 0, 2 };
    plane.normals = { 0,1,0, 0,1,0, 0,1,0, 0,1,0 };
    plane.indices = { 0, 2, 1,  1, 2, 3 };

    bromesh::MeshData ball = bromesh::sphere(1.0f, 8, 6);
    bromesh::translateMesh(ball, 0, 0.5f, 0);
    bromesh::shrinkwrap(ball, plane, bromesh::ShrinkwrapMode::Nearest, 0.0f, 0.0f);

    float maxY = -1e9f, minY = 1e9f;
    for (size_t v = 0; v < ball.vertexCount(); ++v) {
        float y = ball.positions[v*3+1];
        maxY = std::max(maxY, y); minY = std::min(minY, y);
    }
    ASSERT(std::fabs(maxY) < 1e-3f, "sphere verts projected onto y=0 plane (max)");
    ASSERT(std::fabs(minY) < 1e-3f, "sphere verts projected onto y=0 plane (min)");
}

TEST(bake_normal_from_reference) {
    auto low = bromesh::sphere(1.0f, 32, 24);
    bromesh::unwrapUVs(low);
    auto high = bromesh::sphere(1.05f, 64, 48);

    auto nrm = bromesh::bakeNormalsFromReference(low, high, 64, 64, 0.2f);
    ASSERT(nrm.width == 64 && nrm.height == 64, "bake dims");
    ASSERT(nrm.pixels.size() == 64 * 64 * 4, "bake pixel count");

    // Count covered texels and check most encode near +Z in tangent space (b channel ~1)
    int covered = 0, nearUp = 0;
    for (int i = 0; i < 64*64; ++i) {
        if (nrm.pixels[i*4+3] > 0.5f) {
            ++covered;
            if (nrm.pixels[i*4+2] > 0.7f) ++nearUp;
        }
    }
    ASSERT(covered > 100, "bake has coverage");
    ASSERT(nearUp > covered / 2, "most tangent-space normals near +Z");
}

TEST(two_bone_ik_reaches_target) {
    // Skeleton: shoulder at origin, elbow 1 unit along +Y, wrist 1 unit along +Y
    bromesh::Skeleton sk;
    bromesh::Bone s; s.name="s"; s.parent=-1; s.localT[1]=0;
    bromesh::Bone e; e.name="e"; e.parent=0;  e.localT[1]=1;
    bromesh::Bone w; w.name="w"; w.parent=1;  w.localT[1]=1;
    sk.bones = { s, e, w };

    auto pose = bromesh::bindPose(sk);
    float target[3] = { 1.0f, 1.0f, 0.0f };
    bool ok = bromesh::solveTwoBoneIK(sk, pose, 0, 1, 2, target);
    ASSERT(ok, "two-bone IK returns ok");

    // Verify wrist world position is at target
    std::vector<float> world;
    bromesh::computeWorldMatrices(sk, pose, world);
    float wx = world[2*16 + 12], wy = world[2*16 + 13], wz = world[2*16 + 14];
    float dx = wx - target[0], dy = wy - target[1], dz = wz - target[2];
    float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
    ASSERT(dist < 1e-2f, "IK tip reaches target");
}

TEST(fabrik_reaches_target) {
    // 4-bone chain along +Y, each of length 1
    bromesh::Skeleton sk;
    for (int i = 0; i < 4; ++i) {
        bromesh::Bone b;
        b.name = std::string("b") + (char)('0'+i);
        b.parent = (i == 0) ? -1 : (i - 1);
        b.localT[1] = (i == 0) ? 0.0f : 1.0f;
        sk.bones.push_back(b);
    }
    auto pose = bromesh::bindPose(sk);
    std::vector<int> chain = { 0, 1, 2, 3 };
    float target[3] = { 2.0f, 1.0f, 0.0f };
    bool ok = bromesh::solveFABRIK(sk, pose, chain, target, 20, 1e-4f);
    ASSERT(ok, "FABRIK returns ok");

    std::vector<float> world;
    bromesh::computeWorldMatrices(sk, pose, world);
    float wx = world[3*16 + 12], wy = world[3*16 + 13], wz = world[3*16 + 14];
    float dx = wx - target[0], dy = wy - target[1], dz = wz - target[2];
    float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
    ASSERT(dist < 1e-2f, "FABRIK tip reaches target");
}

TEST(retarget_by_name) {
    auto src = makeTwoBoneSkeleton();
    auto dst = makeTwoBoneSkeleton();
    // Change dst order
    std::swap(dst.bones[0], dst.bones[1]);
    // After swap, "child" is at index 0; fix parents
    dst.bones[0].parent = 1; // child has parent root
    dst.bones[1].parent = -1; // root has no parent

    bromesh::Animation anim;
    bromesh::AnimChannel ch;
    ch.boneIndex = 1; // "child" in src
    ch.path = bromesh::AnimChannel::Path::Translation;
    ch.times = { 0.0f }; ch.values = { 1, 2, 3 };
    anim.channels.push_back(ch);
    anim.duration = 0.0f;

    auto remap = bromesh::retargetAnimation(anim, src, dst);
    ASSERT(remap.channels.size() == 1, "retarget preserves channel");
    ASSERT(remap.channels[0].boneIndex == 0, "retargeted to new index of 'child'");
}

TEST(rigify_sockets) {
    bromesh::Skeleton sk;
    for (const char* name : { "DEF-hand.R", "DEF-hand.L", "DEF-head", "DEF-foot.R", "DEF-foot.L" }) {
        bromesh::Bone b; b.name = name; b.parent = -1; sk.bones.push_back(b);
    }
    int added = bromesh::addRigifySockets(sk);
    ASSERT(added >= 5, "rigify helper added at least 5 sockets");
    ASSERT(sk.findSocket("hand.R") >= 0, "hand.R socket present");
    ASSERT(sk.findSocket("head") >= 0, "head socket present");
}

