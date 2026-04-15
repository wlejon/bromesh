#include "test_framework.h"
#include "synthetic_rigs.h"

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

// --- Phase-3: procedural locomotion ---------------------------------------

namespace phase3 {

static bromesh::Skeleton fitCreature(const bromesh::RigSpec& spec,
                                     const bromesh::Landmarks& lm,
                                     const bromesh::MeshData& mesh) {
    return bromesh::fitSkeleton(spec, lm, mesh);
}

static int groundedCount(const std::vector<bromesh::LegChain>& chains) {
    int n = 0;
    for (const auto& c : chains) if (c.grounded) ++n;
    return n;
}

static void sampleTipPositions(const bromesh::Skeleton& sk,
                               const bromesh::Animation& anim,
                               float t,
                               const std::vector<int>& tips,
                               std::vector<float>& outXYZ) {
    auto pose = bromesh::evaluateAnimation(sk, anim, t, /*loop=*/true);
    std::vector<float> world;
    bromesh::computeWorldMatrices(sk, pose, world);
    outXYZ.resize(tips.size() * 3);
    for (size_t i = 0; i < tips.size(); ++i) {
        outXYZ[i*3+0] = world[tips[i]*16 + 12];
        outXYZ[i*3+1] = world[tips[i]*16 + 13];
        outXYZ[i*3+2] = world[tips[i]*16 + 14];
    }
}

} // namespace phase3

TEST(locomotion_identify_legs_humanoid) {
    auto spec = bromesh::builtinHumanoidSpec();
    auto sk = phase3::fitCreature(spec, synthrig::makeHumanoidLandmarks(),
                                  synthrig::makeHumanoidMesh());
    auto chains = bromesh::identifyLegChains(sk, spec);
    ASSERT(chains.size() == 4, "humanoid: 4 IK chains total (2 arms + 2 legs)");
    ASSERT(phase3::groundedCount(chains) == 2, "humanoid: 2 grounded");

    // Grounded ones come first; their tip Z are all ~0 and X are -, + respectively.
    ASSERT(chains[0].grounded && chains[1].grounded, "grounded first");
    ASSERT(chains[0].bones.size() == 3, "leg chain is hip->knee->foot");
    // Each 3-bone chain: parent relations must hold for two-bone IK.
    for (int k = 0; k < 2; ++k) {
        const auto& c = chains[k];
        ASSERT(sk.bones[c.bones[1]].parent == c.bones[0], "mid parent = root");
        ASSERT(sk.bones[c.bones[2]].parent == c.bones[1], "end parent = mid");
    }
}

TEST(locomotion_identify_legs_quadruped) {
    auto spec = bromesh::builtinQuadrupedSpec();
    auto sk = phase3::fitCreature(spec, synthrig::makeQuadrupedLandmarks(),
                                  synthrig::makeQuadrupedMesh());
    auto chains = bromesh::identifyLegChains(sk, spec);
    ASSERT(chains.size() == 4, "quadruped: 4 IK chains (paws only)");
    ASSERT(phase3::groundedCount(chains) == 4, "all 4 grounded");
    // Order: front-L, front-R, hind-L, hind-R (by tip Z desc, then X asc).
    ASSERT(sk.bones[chains[0].tipBone].name == "fl_paw_L", "order[0]=fl_paw_L");
    ASSERT(sk.bones[chains[1].tipBone].name == "fl_paw_R", "order[1]=fl_paw_R");
    ASSERT(sk.bones[chains[2].tipBone].name == "hl_paw_L", "order[2]=hl_paw_L");
    ASSERT(sk.bones[chains[3].tipBone].name == "hl_paw_R", "order[3]=hl_paw_R");
}

TEST(locomotion_identify_legs_hexapod) {
    auto spec = bromesh::builtinHexapodSpec();
    auto sk = phase3::fitCreature(spec, synthrig::makeHexapodLandmarks(),
                                  synthrig::makeHexapodMesh());
    auto chains = bromesh::identifyLegChains(sk, spec);
    ASSERT(chains.size() == 6, "hexapod: 6 IK chains");
    ASSERT(phase3::groundedCount(chains) == 6, "all 6 grounded");
    const char* expected[6] = {
        "front_tibia_L", "front_tibia_R",
        "mid_tibia_L",   "mid_tibia_R",
        "rear_tibia_L",  "rear_tibia_R",
    };
    for (int i = 0; i < 6; ++i) {
        ASSERT(sk.bones[chains[i].tipBone].name == expected[i], "hexapod leg order");
    }
}

TEST(locomotion_identify_legs_octopod) {
    auto spec = bromesh::builtinOctopodSpec();
    auto sk = phase3::fitCreature(spec, synthrig::makeOctopodLandmarks(),
                                  synthrig::makeOctopodMesh());
    auto chains = bromesh::identifyLegChains(sk, spec);
    ASSERT(chains.size() == 8, "octopod: 8 IK chains");
    ASSERT(phase3::groundedCount(chains) == 8, "all 8 grounded");
    // Each octopod chain is only 2 bones (upper, lower) → FABRIK path at runtime.
    for (const auto& c : chains) {
        ASSERT(c.bones.size() == 2, "octopod arm chain has 2 bones");
    }
}

TEST(locomotion_default_gait_sizes) {
    ASSERT(bromesh::defaultGait(2).phases.size() == 2, "biped gait size");
    ASSERT(bromesh::defaultGait(4).phases.size() == 4, "quadruped gait size");
    ASSERT(bromesh::defaultGait(6).phases.size() == 6, "hexapod gait size");
    ASSERT(bromesh::defaultGait(8).phases.size() == 8, "octopod gait size");
    ASSERT(bromesh::defaultGait(3).phases.empty(), "odd leg count has no default");
    ASSERT(bromesh::defaultGait(2).dutyFactor > 0.0f &&
           bromesh::defaultGait(2).dutyFactor < 1.0f, "duty in (0,1)");
}

TEST(locomotion_cycle_humanoid_closes) {
    auto spec = bromesh::builtinHumanoidSpec();
    auto sk = phase3::fitCreature(spec, synthrig::makeHumanoidLandmarks(),
                                  synthrig::makeHumanoidMesh());
    bromesh::LocomotionParams p;
    p.cycleDuration = 1.0f;
    p.strideLength = 0.25f;
    auto anim = bromesh::generateLocomotionCycle(sk, spec, p);
    ASSERT(std::fabs(anim.duration - 1.0f) < 1e-5f, "duration == cycleDuration");
    ASSERT(!anim.channels.empty(), "channels emitted");

    // Loop closure: tip positions at t=0 and t=duration must match.
    auto chains = bromesh::identifyLegChains(sk, spec);
    std::vector<int> tips;
    for (const auto& c : chains) if (c.grounded) tips.push_back(c.tipBone);

    std::vector<float> p0, p1;
    phase3::sampleTipPositions(sk, anim, 0.0f,            tips, p0);
    phase3::sampleTipPositions(sk, anim, anim.duration,   tips, p1);
    float maxD = 0.0f;
    for (size_t i = 0; i < p0.size(); ++i) {
        float d = std::fabs(p0[i] - p1[i]);
        if (d > maxD) maxD = d;
    }
    ASSERT(maxD < 1e-3f, "loop closure: t=0 matches t=duration");
}

TEST(locomotion_feet_actually_move) {
    auto spec = bromesh::builtinHumanoidSpec();
    auto sk = phase3::fitCreature(spec, synthrig::makeHumanoidLandmarks(),
                                  synthrig::makeHumanoidMesh());
    bromesh::LocomotionParams p;
    p.cycleDuration = 1.0f;
    p.strideLength = 0.25f;
    auto anim = bromesh::generateLocomotionCycle(sk, spec, p);

    auto chains = bromesh::identifyLegChains(sk, spec);
    std::vector<int> tips;
    for (const auto& c : chains) if (c.grounded) tips.push_back(c.tipBone);

    std::vector<float> p0, pHalf;
    phase3::sampleTipPositions(sk, anim, 0.0f, tips, p0);
    phase3::sampleTipPositions(sk, anim, 0.5f, tips, pHalf);

    // At t=0.5 (half cycle), with biped [0, 0.5] phases, one foot is planted
    // near the rear and one is near the front: both should have moved from
    // their bind positions by a visible fraction of stride.
    float maxD = 0.0f;
    for (size_t i = 0; i < p0.size(); ++i) {
        float d = std::fabs(p0[i] - pHalf[i]);
        if (d > maxD) maxD = d;
    }
    ASSERT(maxD > p.strideLength * 0.15f, "feet move across the cycle");
}

TEST(locomotion_deterministic) {
    auto spec = bromesh::builtinQuadrupedSpec();
    auto sk = phase3::fitCreature(spec, synthrig::makeQuadrupedLandmarks(),
                                  synthrig::makeQuadrupedMesh());
    bromesh::LocomotionParams p;
    p.cycleDuration = 0.8f;
    auto a = bromesh::generateLocomotionCycle(sk, spec, p);
    auto b = bromesh::generateLocomotionCycle(sk, spec, p);
    ASSERT(a.channels.size() == b.channels.size(), "same channel count");
    bool allEqual = true;
    for (size_t i = 0; i < a.channels.size(); ++i) {
        if (a.channels[i].values != b.channels[i].values) { allEqual = false; break; }
        if (a.channels[i].times  != b.channels[i].times)  { allEqual = false; break; }
    }
    ASSERT(allEqual, "channel data byte-identical across calls");
}

TEST(locomotion_cycle_all_creatures_smoke) {
    struct Case {
        const char* name;
        bromesh::RigSpec spec;
        bromesh::Skeleton sk;
    };
    std::vector<Case> cases;
    {
        Case c; c.name = "humanoid"; c.spec = bromesh::builtinHumanoidSpec();
        c.sk = phase3::fitCreature(c.spec, synthrig::makeHumanoidLandmarks(),
                                   synthrig::makeHumanoidMesh());
        cases.push_back(std::move(c));
    }
    {
        Case c; c.name = "quadruped"; c.spec = bromesh::builtinQuadrupedSpec();
        c.sk = phase3::fitCreature(c.spec, synthrig::makeQuadrupedLandmarks(),
                                   synthrig::makeQuadrupedMesh());
        cases.push_back(std::move(c));
    }
    {
        Case c; c.name = "hexapod"; c.spec = bromesh::builtinHexapodSpec();
        c.sk = phase3::fitCreature(c.spec, synthrig::makeHexapodLandmarks(),
                                   synthrig::makeHexapodMesh());
        cases.push_back(std::move(c));
    }
    {
        Case c; c.name = "octopod"; c.spec = bromesh::builtinOctopodSpec();
        c.sk = phase3::fitCreature(c.spec, synthrig::makeOctopodLandmarks(),
                                   synthrig::makeOctopodMesh());
        cases.push_back(std::move(c));
    }

    bromesh::LocomotionParams p;
    p.cycleDuration = 1.0f;
    p.strideLength = 0.20f;
    p.footLiftHeight = 0.06f;
    p.keyframesPerCycle = 16;

    for (const auto& c : cases) {
        auto anim = bromesh::generateLocomotionCycle(c.sk, c.spec, p);
        ASSERT(anim.duration > 0.0f, "anim has duration");
        ASSERT(!anim.channels.empty(), "anim has channels");

        // Sample at multiple times and verify pose values finite, and
        // skinning matrix composition succeeds without NaN.
        bool allFinite = true;
        for (int i = 0; i <= 8; ++i) {
            float t = (float)i / 8.0f * anim.duration;
            auto pose = bromesh::evaluateAnimation(c.sk, anim, t, true);
            for (float v : pose.data) {
                if (!(v == v)) { allFinite = false; break; }
            }
            if (!allFinite) break;
            std::vector<float> skm;
            bromesh::computeSkinningMatrices(c.sk, pose, skm);
            for (float v : skm) {
                if (!(v == v)) { allFinite = false; break; }
            }
            if (!allFinite) break;
        }
        ASSERT(allFinite, "locomotion pose values finite across cycle");
    }
}

