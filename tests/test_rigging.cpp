#include "test_framework.h"

// --- Auto-rigging tests ----------------------------------------------------

static bromesh::MeshData makeSyntheticHumanoid() {
    // Stacked boxes forming a stick-figure humanoid.
    // Coordinate convention: +Y up, +Z forward, +X right (so "Left" body
    // side is at negative X — we pick a convention and stick to it below
    // when placing landmarks).
    auto boxAt = [](float cx, float cy, float cz,
                    float hx, float hy, float hz) {
        auto m = bromesh::box(hx, hy, hz);
        bromesh::translateMesh(m, cx, cy, cz);
        return m;
    };
    std::vector<bromesh::MeshData> parts;
    // Torso
    parts.push_back(boxAt(0.0f,  0.25f, 0.0f, 0.18f, 0.30f, 0.10f));
    // Head
    parts.push_back(boxAt(0.0f,  0.70f, 0.0f, 0.12f, 0.12f, 0.12f));
    // Left arm (negative X)
    parts.push_back(boxAt(-0.35f, 0.45f, 0.0f, 0.18f, 0.06f, 0.06f));
    // Right arm (positive X)
    parts.push_back(boxAt( 0.35f, 0.45f, 0.0f, 0.18f, 0.06f, 0.06f));
    // Left leg
    parts.push_back(boxAt(-0.09f,-0.45f, 0.0f, 0.06f, 0.45f, 0.06f));
    // Right leg
    parts.push_back(boxAt( 0.09f,-0.45f, 0.0f, 0.06f, 0.45f, 0.06f));
    return bromesh::mergeMeshes(parts);
}

static bromesh::Landmarks makeHumanoidLandmarks() {
    bromesh::Landmarks lm;
    // Axial
    lm.set("pelvis",    0.00f, -0.05f, 0.00f);
    lm.set("chest",     0.00f,  0.45f, 0.00f);
    lm.set("neck_base", 0.00f,  0.55f, 0.00f);
    lm.set("crown",     0.00f,  0.80f, 0.00f);
    // Arms (L at -X, R at +X; names follow spec mirror pairs)
    lm.set("shoulder_L",-0.18f, 0.50f, 0.00f);
    lm.set("shoulder_R", 0.18f, 0.50f, 0.00f);
    lm.set("elbow_L",   -0.35f, 0.45f, 0.00f);
    lm.set("elbow_R",    0.35f, 0.45f, 0.00f);
    lm.set("wrist_L",   -0.50f, 0.45f, 0.00f);
    lm.set("wrist_R",    0.50f, 0.45f, 0.00f);
    // Legs
    lm.set("hip_L",     -0.09f, -0.05f, 0.00f);
    lm.set("hip_R",      0.09f, -0.05f, 0.00f);
    lm.set("knee_L",    -0.09f, -0.45f, 0.00f);
    lm.set("knee_R",     0.09f, -0.45f, 0.00f);
    lm.set("ankle_L",   -0.09f, -0.88f, 0.00f);
    lm.set("ankle_R",    0.09f, -0.88f, 0.00f);
    lm.set("toe_L",     -0.09f, -0.88f, 0.05f);
    lm.set("toe_R",      0.09f, -0.88f, 0.05f);
    return lm;
}

TEST(rig_spec_humanoid_shape) {
    auto spec = bromesh::builtinHumanoidSpec();
    ASSERT(spec.name == "humanoid", "spec name");
    ASSERT(spec.symmetric, "humanoid is symmetric");
    ASSERT(spec.bones.size() == 22, "humanoid has 22 bones");
    ASSERT(spec.landmarks.size() == 18, "humanoid has 18 landmark decls");
    ASSERT(spec.sockets.size() == 3, "humanoid has 3 default sockets");

    // Every non-root bone's parent must be present by name.
    std::unordered_set<std::string> names;
    for (const auto& b : spec.bones) names.insert(b.name);
    for (const auto& b : spec.bones) {
        if (b.parent.empty()) continue;
        ASSERT(names.count(b.parent) == 1, "parent name resolves");
    }
}

TEST(rig_spec_json_roundtrip) {
    auto spec = bromesh::builtinHumanoidSpec();
    std::string js = bromesh::serializeRigSpecJSON(spec);
    // With tinygltf-backed json.hpp the string is non-empty; without it we
    // fall through — accept either to keep the test robust to the build
    // configuration.
    if (js.empty()) return;
    auto parsed = bromesh::parseRigSpecJSON(js);
    ASSERT(parsed.name == spec.name, "json roundtrip name");
    ASSERT(parsed.bones.size() == spec.bones.size(), "json roundtrip bone count");
    ASSERT(parsed.landmarks.size() == spec.landmarks.size(), "json roundtrip landmark count");
    ASSERT(parsed.sockets.size() == spec.sockets.size(), "json roundtrip socket count");
}

TEST(skeleton_fit_humanoid) {
    auto spec = bromesh::builtinHumanoidSpec();
    auto lm = makeHumanoidLandmarks();
    auto missing = bromesh::missingLandmarks(spec, lm);
    ASSERT(missing.empty(), "no missing landmarks");

    auto mesh = makeSyntheticHumanoid();
    auto skel = bromesh::fitSkeleton(spec, lm, mesh);
    ASSERT(skel.bones.size() == spec.bones.size(), "bone count matches spec");
    ASSERT(skel.sockets.size() == 3, "sockets emitted");

    // Parents topologically precede children.
    for (size_t i = 0; i < skel.bones.size(); ++i) {
        ASSERT(skel.bones[i].parent < (int)i, "parent precedes child");
    }

    // At least one root bone.
    bool hasRoot = false;
    for (const auto& b : skel.bones) if (b.parent == -1) { hasRoot = true; break; }
    ASSERT(hasRoot, "skeleton has a root");
}

TEST(auto_rig_end_to_end) {
    auto spec = bromesh::builtinHumanoidSpec();
    auto lm = makeHumanoidLandmarks();
    auto mesh = makeSyntheticHumanoid();

    bromesh::VoxelBindOptions opts;
    opts.maxResolution = 48; // keep the test fast
    auto r = bromesh::autoRig(mesh, spec, lm, opts);

    ASSERT(r.missingLandmarks.empty(), "no missing landmarks");
    ASSERT(r.skeleton.bones.size() == spec.bones.size(), "skeleton bone count");
    ASSERT(r.skin.boneCount == r.skeleton.bones.size(), "skin boneCount");
    ASSERT(r.skin.boneWeights.size() == mesh.vertexCount() * 4, "weights sized");
    ASSERT(r.skin.boneIndices.size() == mesh.vertexCount() * 4, "indices sized");
    ASSERT(r.skin.inverseBindMatrices.size() == r.skeleton.bones.size() * 16,
           "inverseBind sized");

    // Every vertex should have weights summing to ~1 and at least one influence.
    size_t bad = 0, orphan = 0;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float sum = 0.0f; int nz = 0;
        for (int k = 0; k < 4; ++k) {
            float w = r.skin.boneWeights[v * 4 + k];
            if (w != w) { bad++; break; }
            sum += w;
            if (w > 0.0f) ++nz;
        }
        if (std::fabs(sum - 1.0f) > 1e-3f) ++bad;
        if (nz == 0) ++orphan;
    }
    ASSERT(bad == 0, "all vertices have valid weight sum");
    ASSERT(orphan == 0, "no orphan vertices");

    // Bind-pose skinning must leave positions unchanged. applySkinning
    // expects *world* matrices per bone (and multiplies by inverseBind
    // internally via SkinData), so we pass computeWorldMatrices's output
    // directly rather than computeSkinningMatrices's.
    auto pose = bromesh::bindPose(r.skeleton);
    std::vector<float> world;
    bromesh::computeWorldMatrices(r.skeleton, pose, world);
    auto beforeMesh = mesh;
    bromesh::applySkinning(beforeMesh, r.skin, world.data());
    float maxDelta = 0.0f;
    for (size_t i = 0; i < mesh.positions.size(); ++i) {
        float d = std::fabs(beforeMesh.positions[i] - mesh.positions[i]);
        if (d > maxDelta) maxDelta = d;
    }
    ASSERT(maxDelta < 1e-3f, "bind pose skinning is identity");
}

TEST(auto_rig_deterministic) {
    auto spec = bromesh::builtinHumanoidSpec();
    auto lm = makeHumanoidLandmarks();
    auto mesh = makeSyntheticHumanoid();
    bromesh::VoxelBindOptions opts; opts.maxResolution = 48;
    auto a = bromesh::autoRig(mesh, spec, lm, opts);
    auto b = bromesh::autoRig(mesh, spec, lm, opts);
    ASSERT(a.skin.boneWeights == b.skin.boneWeights, "weights deterministic");
    ASSERT(a.skin.boneIndices == b.skin.boneIndices, "indices deterministic");
}

// --- Phase-2: non-humanoid rig specs ---------------------------------------
//
// The tests below verify that the data-driven RigSpec schema generalizes to
// quadruped / hexapod / octopod rigs without any change to fitSkeleton,
// voxelBindWeights, or autoRig. Each creature gets a synthetic box-assembly
// mesh, a full landmark dict, and the same four assertions: spec shape, JSON
// roundtrip, skeleton fit, and end-to-end autoRig weight validity.

namespace phase2 {

static bromesh::MeshData boxAt(float cx, float cy, float cz,
                               float hx, float hy, float hz) {
    auto m = bromesh::box(hx, hy, hz);
    bromesh::translateMesh(m, cx, cy, cz);
    return m;
}

// ---- Quadruped ------------------------------------------------------------

static bromesh::MeshData makeSyntheticQuadruped() {
    std::vector<bromesh::MeshData> parts;
    // Torso (body runs along Z).
    parts.push_back(boxAt(0.00f, 0.45f, 0.00f, 0.12f, 0.12f, 0.30f));
    // Head.
    parts.push_back(boxAt(0.00f, 0.55f, 0.40f, 0.08f, 0.08f, 0.10f));
    // Tail.
    parts.push_back(boxAt(0.00f, 0.50f,-0.45f, 0.03f, 0.03f, 0.15f));
    // Front legs (z = +0.20).
    parts.push_back(boxAt(-0.08f, 0.18f, 0.20f, 0.04f, 0.27f, 0.04f));
    parts.push_back(boxAt( 0.08f, 0.18f, 0.20f, 0.04f, 0.27f, 0.04f));
    // Hind legs (z = -0.20).
    parts.push_back(boxAt(-0.08f, 0.18f,-0.20f, 0.04f, 0.27f, 0.04f));
    parts.push_back(boxAt( 0.08f, 0.18f,-0.20f, 0.04f, 0.27f, 0.04f));
    return bromesh::mergeMeshes(parts);
}

static bromesh::Landmarks makeQuadrupedLandmarks() {
    bromesh::Landmarks lm;
    lm.set("pelvis",     0.00f, 0.45f, -0.22f);
    lm.set("chest",      0.00f, 0.45f,  0.22f);
    lm.set("neck_base",  0.00f, 0.52f,  0.28f);
    lm.set("crown",      0.00f, 0.58f,  0.40f);
    lm.set("muzzle",     0.00f, 0.55f,  0.50f);
    lm.set("tail_base",  0.00f, 0.50f, -0.30f);
    lm.set("tail_tip",   0.00f, 0.50f, -0.60f);
    lm.set("fshoulder_L",-0.10f, 0.42f,  0.22f);
    lm.set("fshoulder_R", 0.10f, 0.42f,  0.22f);
    lm.set("felbow_L",   -0.08f, 0.25f,  0.20f);
    lm.set("felbow_R",    0.08f, 0.25f,  0.20f);
    lm.set("fpaw_L",     -0.08f,-0.08f,  0.20f);
    lm.set("fpaw_R",      0.08f,-0.08f,  0.20f);
    lm.set("hip_L",      -0.10f, 0.42f, -0.20f);
    lm.set("hip_R",       0.10f, 0.42f, -0.20f);
    lm.set("hknee_L",    -0.08f, 0.25f, -0.20f);
    lm.set("hknee_R",     0.08f, 0.25f, -0.20f);
    lm.set("hpaw_L",     -0.08f,-0.08f, -0.20f);
    lm.set("hpaw_R",      0.08f,-0.08f, -0.20f);
    return lm;
}

// ---- Hexapod --------------------------------------------------------------

static bromesh::MeshData makeSyntheticHexapod() {
    std::vector<bromesh::MeshData> parts;
    // Abdomen (rear segment), thorax (mid), head (front).
    parts.push_back(boxAt(0.00f, 0.20f,-0.20f, 0.08f, 0.08f, 0.15f));
    parts.push_back(boxAt(0.00f, 0.20f, 0.00f, 0.10f, 0.10f, 0.12f));
    parts.push_back(boxAt(0.00f, 0.20f, 0.18f, 0.06f, 0.06f, 0.06f));
    // Six legs attached to the thorax, three per side, roughly at z = +0.08 / 0 / -0.08.
    for (int side = 0; side < 2; ++side) {
        float sx = (side == 0) ? -1.0f : 1.0f;
        for (int row = 0; row < 3; ++row) {
            float z = 0.08f - row * 0.08f;
            parts.push_back(boxAt(sx * 0.12f, 0.10f, z, 0.03f, 0.15f, 0.03f));
        }
    }
    return bromesh::mergeMeshes(parts);
}

static bromesh::Landmarks makeHexapodLandmarks() {
    bromesh::Landmarks lm;
    lm.set("abdomen_tip", 0.00f, 0.20f, -0.35f);
    lm.set("abdomen",     0.00f, 0.20f, -0.18f);
    lm.set("thorax",      0.00f, 0.20f,  0.05f);
    lm.set("head",        0.00f, 0.20f,  0.15f);
    lm.set("head_tip",    0.00f, 0.20f,  0.25f);
    const char* prefixes[3] = { "front", "mid", "rear" };
    float rowZ[3] = { 0.08f, 0.00f, -0.08f };
    for (int i = 0; i < 3; ++i) {
        float z = rowZ[i];
        std::string p = prefixes[i];
        lm.set((p + "_hip_L").c_str(),  -0.10f, 0.18f, z);
        lm.set((p + "_hip_R").c_str(),   0.10f, 0.18f, z);
        lm.set((p + "_knee_L").c_str(), -0.18f, 0.12f, z);
        lm.set((p + "_knee_R").c_str(),  0.18f, 0.12f, z);
        lm.set((p + "_foot_L").c_str(), -0.22f, 0.00f, z);
        lm.set((p + "_foot_R").c_str(),  0.22f, 0.00f, z);
    }
    return lm;
}

// ---- Octopod --------------------------------------------------------------

static bromesh::MeshData makeSyntheticOctopod() {
    std::vector<bromesh::MeshData> parts;
    // Round-ish body and a small "head" bump on top.
    parts.push_back(boxAt(0.00f, 0.25f, 0.00f, 0.14f, 0.08f, 0.14f));
    parts.push_back(boxAt(0.00f, 0.35f, 0.00f, 0.06f, 0.04f, 0.06f));
    // Eight arms splay radially from the body; approximate as boxes at
    // arm_a (z=+0.15), arm_b (z=+0.05), arm_c (z=-0.05), arm_d (z=-0.15),
    // each with L/R pair on X axis.
    float armZ[4] = { 0.15f, 0.05f, -0.05f, -0.15f };
    for (int i = 0; i < 4; ++i) {
        float z = armZ[i];
        parts.push_back(boxAt(-0.22f, 0.15f, z, 0.08f, 0.03f, 0.03f));
        parts.push_back(boxAt( 0.22f, 0.15f, z, 0.08f, 0.03f, 0.03f));
    }
    return bromesh::mergeMeshes(parts);
}

static bromesh::Landmarks makeOctopodLandmarks() {
    bromesh::Landmarks lm;
    lm.set("body", 0.00f, 0.25f, 0.00f);
    lm.set("head", 0.00f, 0.38f, 0.00f);
    const char* arms[4] = { "arm_a", "arm_b", "arm_c", "arm_d" };
    float armZ[4] = { 0.15f, 0.05f, -0.05f, -0.15f };
    for (int i = 0; i < 4; ++i) {
        float z = armZ[i];
        std::string a = arms[i];
        lm.set((a + "_hip_L").c_str(), -0.14f, 0.22f, z);
        lm.set((a + "_hip_R").c_str(),  0.14f, 0.22f, z);
        lm.set((a + "_mid_L").c_str(), -0.22f, 0.15f, z);
        lm.set((a + "_mid_R").c_str(),  0.22f, 0.15f, z);
        lm.set((a + "_tip_L").c_str(), -0.30f, 0.10f, z);
        lm.set((a + "_tip_R").c_str(),  0.30f, 0.10f, z);
    }
    return lm;
}

// ---- Parameterized assertions --------------------------------------------

static void assertParentsResolve(const bromesh::RigSpec& spec) {
    std::unordered_set<std::string> names;
    for (const auto& b : spec.bones) names.insert(b.name);
    for (const auto& b : spec.bones) {
        if (b.parent.empty()) continue;
        ASSERT(names.count(b.parent) == 1, "parent name resolves");
    }
}

static void assertJsonRoundtrip(const bromesh::RigSpec& spec) {
    std::string js = bromesh::serializeRigSpecJSON(spec);
    if (js.empty()) return; // JSON support not compiled in
    auto parsed = bromesh::parseRigSpecJSON(js);
    ASSERT(parsed.name == spec.name, "json roundtrip name");
    ASSERT(parsed.bones.size() == spec.bones.size(), "json roundtrip bone count");
    ASSERT(parsed.landmarks.size() == spec.landmarks.size(), "json roundtrip landmark count");
    ASSERT(parsed.sockets.size() == spec.sockets.size(), "json roundtrip socket count");
}

static void assertEndToEnd(const bromesh::RigSpec& spec,
                           const bromesh::Landmarks& lm,
                           const bromesh::MeshData& mesh) {
    auto missing = bromesh::missingLandmarks(spec, lm);
    ASSERT(missing.empty(), "no missing landmarks");

    bromesh::VoxelBindOptions opts;
    opts.maxResolution = 48;
    auto r = bromesh::autoRig(mesh, spec, lm, opts);

    ASSERT(r.missingLandmarks.empty(), "autoRig: no missing landmarks");
    ASSERT(r.skeleton.bones.size() == spec.bones.size(), "skeleton bone count");
    ASSERT(r.skin.boneCount == r.skeleton.bones.size(), "skin boneCount");
    ASSERT(r.skin.boneWeights.size() == mesh.vertexCount() * 4, "weights sized");

    // Parents topologically precede children.
    for (size_t i = 0; i < r.skeleton.bones.size(); ++i) {
        ASSERT(r.skeleton.bones[i].parent < (int)i, "parent precedes child");
    }

    // Weight sanity: sum ~1, no NaNs, no orphan vertices.
    size_t bad = 0, orphan = 0;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float sum = 0.0f; int nz = 0;
        for (int k = 0; k < 4; ++k) {
            float w = r.skin.boneWeights[v * 4 + k];
            if (w != w) { bad++; break; }
            sum += w;
            if (w > 0.0f) ++nz;
        }
        if (std::fabs(sum - 1.0f) > 1e-3f) ++bad;
        if (nz == 0) ++orphan;
    }
    ASSERT(bad == 0, "all vertices have valid weight sum");
    ASSERT(orphan == 0, "no orphan vertices");

    // Bind-pose skinning is identity.
    auto pose = bromesh::bindPose(r.skeleton);
    std::vector<float> world;
    bromesh::computeWorldMatrices(r.skeleton, pose, world);
    auto skinned = mesh;
    bromesh::applySkinning(skinned, r.skin, world.data());
    float maxDelta = 0.0f;
    for (size_t i = 0; i < mesh.positions.size(); ++i) {
        float d = std::fabs(skinned.positions[i] - mesh.positions[i]);
        if (d > maxDelta) maxDelta = d;
    }
    ASSERT(maxDelta < 1e-3f, "bind pose skinning is identity");
}

} // namespace phase2

TEST(rig_spec_quadruped_shape) {
    auto spec = bromesh::builtinQuadrupedSpec();
    ASSERT(spec.name == "quadruped", "spec name");
    ASSERT(spec.symmetric, "quadruped is symmetric");
    ASSERT(spec.bones.size() == 20, "quadruped has 20 bones");
    ASSERT(spec.landmarks.size() == 19, "quadruped has 19 landmark decls");
    ASSERT(spec.sockets.size() == 2, "quadruped has 2 default sockets");
    phase2::assertParentsResolve(spec);
}

TEST(rig_spec_quadruped_json_roundtrip) {
    phase2::assertJsonRoundtrip(bromesh::builtinQuadrupedSpec());
}

TEST(auto_rig_quadruped_end_to_end) {
    phase2::assertEndToEnd(bromesh::builtinQuadrupedSpec(),
                           phase2::makeQuadrupedLandmarks(),
                           phase2::makeSyntheticQuadruped());
}

TEST(rig_spec_hexapod_shape) {
    auto spec = bromesh::builtinHexapodSpec();
    ASSERT(spec.name == "hexapod", "spec name");
    ASSERT(spec.symmetric, "hexapod is symmetric");
    ASSERT(spec.bones.size() == 21, "hexapod has 21 bones");
    ASSERT(spec.landmarks.size() == 23, "hexapod has 23 landmark decls");
    ASSERT(spec.sockets.size() == 1, "hexapod has 1 default socket");
    phase2::assertParentsResolve(spec);
}

TEST(rig_spec_hexapod_json_roundtrip) {
    phase2::assertJsonRoundtrip(bromesh::builtinHexapodSpec());
}

TEST(auto_rig_hexapod_end_to_end) {
    phase2::assertEndToEnd(bromesh::builtinHexapodSpec(),
                           phase2::makeHexapodLandmarks(),
                           phase2::makeSyntheticHexapod());
}

TEST(rig_spec_octopod_shape) {
    auto spec = bromesh::builtinOctopodSpec();
    ASSERT(spec.name == "octopod", "spec name");
    ASSERT(spec.symmetric, "octopod is symmetric");
    ASSERT(spec.bones.size() == 18, "octopod has 18 bones");
    ASSERT(spec.landmarks.size() == 26, "octopod has 26 landmark decls");
    ASSERT(spec.sockets.size() == 1, "octopod has 1 default socket");
    phase2::assertParentsResolve(spec);
}

TEST(rig_spec_octopod_json_roundtrip) {
    phase2::assertJsonRoundtrip(bromesh::builtinOctopodSpec());
}

TEST(auto_rig_octopod_end_to_end) {
    phase2::assertEndToEnd(bromesh::builtinOctopodSpec(),
                           phase2::makeOctopodLandmarks(),
                           phase2::makeSyntheticOctopod());
}

// --- Phase-5: heuristic quadruped landmark detection ---------------------

TEST(detect_landmarks_quadruped_completeness) {
    auto mesh = phase2::makeSyntheticQuadruped();
    auto lm = bromesh::detectQuadrupedLandmarks(mesh);
    auto spec = bromesh::builtinQuadrupedSpec();
    auto missing = bromesh::missingLandmarks(spec, lm);
    ASSERT(missing.empty(), "all 19 quadruped landmarks detected");
}

TEST(detect_landmarks_quadruped_near_reference) {
    auto mesh = phase2::makeSyntheticQuadruped();
    auto detected = bromesh::detectQuadrupedLandmarks(mesh);
    auto reference = phase2::makeQuadrupedLandmarks();

    // Quadruped body is elongated along forward; tolerance tracks the
    // largest extent so the bar doesn't vary wildly with species shape.
    auto bbox = bromesh::computeBBox(mesh);
    float scale = std::max({bbox.extentX(), bbox.extentY(), bbox.extentZ()});
    float tol = 0.15f * scale;

    int checked = 0;
    for (const auto& [name, ref] : reference.points) {
        if (!detected.has(name)) continue;
        auto d = detected.points[name];
        float dx = d[0]-ref[0], dy = d[1]-ref[1], dz = d[2]-ref[2];
        float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
        ASSERT(dist < tol, name.c_str());
        ++checked;
    }
    ASSERT(checked == 19, "checked all 19 landmarks");
}

TEST(detect_landmarks_quadruped_end_to_end) {
    auto mesh = phase2::makeSyntheticQuadruped();
    auto spec = bromesh::builtinQuadrupedSpec();
    auto lm = bromesh::detectQuadrupedLandmarks(mesh);

    bromesh::VoxelBindOptions opts; opts.maxResolution = 48;
    auto r = bromesh::autoRig(mesh, spec, lm, opts);
    ASSERT(r.missingLandmarks.empty(), "no missing landmarks");
    ASSERT(r.skeleton.bones.size() == spec.bones.size(), "bone count");

    auto pose = bromesh::bindPose(r.skeleton);
    std::vector<float> world;
    bromesh::computeWorldMatrices(r.skeleton, pose, world);
    auto skinned = mesh;
    bromesh::applySkinning(skinned, r.skin, world.data());
    float maxDelta = 0.0f;
    for (size_t i = 0; i < mesh.positions.size(); ++i) {
        float d = std::fabs(skinned.positions[i] - mesh.positions[i]);
        if (d > maxDelta) maxDelta = d;
    }
    ASSERT(maxDelta < 1e-3f, "bind pose skinning is identity");
}

TEST(detect_landmarks_quadruped_deterministic) {
    auto mesh = phase2::makeSyntheticQuadruped();
    auto a = bromesh::detectQuadrupedLandmarks(mesh);
    auto b = bromesh::detectQuadrupedLandmarks(mesh);
    ASSERT(a.points.size() == b.points.size(), "same count");
    for (const auto& [name, pa] : a.points) {
        ASSERT(b.has(name), name.c_str());
        auto pb = b.points[name];
        ASSERT(pa[0] == pb[0] && pa[1] == pb[1] && pa[2] == pb[2],
               "same position both calls");
    }
}

TEST(auto_rig_quadruped_deterministic) {
    auto spec = bromesh::builtinQuadrupedSpec();
    auto lm = phase2::makeQuadrupedLandmarks();
    auto mesh = phase2::makeSyntheticQuadruped();
    bromesh::VoxelBindOptions opts; opts.maxResolution = 48;
    auto a = bromesh::autoRig(mesh, spec, lm, opts);
    auto b = bromesh::autoRig(mesh, spec, lm, opts);
    ASSERT(a.skin.boneWeights == b.skin.boneWeights, "weights deterministic");
    ASSERT(a.skin.boneIndices == b.skin.boneIndices, "indices deterministic");
}

// --- Phase-4: weight post-processing (smoothing + outlier rejection) ----

// Helper: compute mean neighbor-to-neighbor weight difference across the
// mesh. Lower = smoother (voxel-grid stairsteps show up as high values).
static float meshWeightRoughness(const bromesh::MeshData& mesh,
                                 const bromesh::SkinData& skin) {
    const size_t nV = mesh.vertexCount();
    if (nV == 0) return 0.0f;
    const size_t K = 4;
    // Build adjacency.
    std::vector<std::unordered_set<uint32_t>> adj(nV);
    for (size_t t = 0; t + 2 < mesh.indices.size(); t += 3) {
        uint32_t a = mesh.indices[t], b = mesh.indices[t+1], c = mesh.indices[t+2];
        if (a>=nV||b>=nV||c>=nV) continue;
        adj[a].insert(b); adj[a].insert(c);
        adj[b].insert(a); adj[b].insert(c);
        adj[c].insert(a); adj[c].insert(b);
    }
    // Dense per-vertex weight table.
    const size_t B = skin.boneCount;
    std::vector<float> dense(nV * B, 0.0f);
    for (size_t v = 0; v < nV; ++v) {
        for (size_t k = 0; k < K; ++k) {
            float w = skin.boneWeights[v*K+k];
            uint32_t bi = skin.boneIndices[v*K+k];
            if (w > 0.0f && bi < B) dense[v*B + bi] += w;
        }
    }
    double total = 0.0; size_t edges = 0;
    for (size_t v = 0; v < nV; ++v) {
        for (uint32_t n : adj[v]) {
            if (n <= v) continue;
            for (size_t bi = 0; bi < B; ++bi) {
                total += std::fabs(dense[v*B + bi] - dense[(size_t)n*B + bi]);
            }
            ++edges;
        }
    }
    return edges ? float(total / (double)edges) : 0.0f;
}

TEST(post_process_preserves_sum_to_one) {
    auto spec = bromesh::builtinHumanoidSpec();
    auto lm = makeHumanoidLandmarks();
    auto mesh = makeSyntheticHumanoid();
    bromesh::VoxelBindOptions opts; opts.maxResolution = 48;
    auto r = bromesh::autoRig(mesh, spec, lm, opts);
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float sum = 0.0f;
        for (int k = 0; k < 4; ++k) sum += r.skin.boneWeights[v*4+k];
        ASSERT(std::fabs(sum - 1.0f) < 1e-3f, "weights sum to 1 after post-process");
    }
}

TEST(post_process_reduces_roughness) {
    auto spec = bromesh::builtinHumanoidSpec();
    auto lm = makeHumanoidLandmarks();
    auto mesh = makeSyntheticHumanoid();

    bromesh::VoxelBindOptions raw;  raw.maxResolution = 48; raw.smoothIterations = 0;
    bromesh::VoxelBindOptions smth; smth.maxResolution = 48; smth.smoothIterations = 6;
    auto rRaw = bromesh::autoRig(mesh, spec, lm, raw);
    auto rSmth = bromesh::autoRig(mesh, spec, lm, smth);

    float rough0 = meshWeightRoughness(mesh, rRaw.skin);
    float rough1 = meshWeightRoughness(mesh, rSmth.skin);
    ASSERT(rough1 <= rough0, "smoothing does not increase roughness");
    // Expect a meaningful reduction — at least 15%. Voxel-bind on this mesh
    // produces visible stepping so smoothing should win clearly.
    ASSERT(rough1 < rough0 * 0.85f, "smoothing reduces roughness by >=15%");
}

TEST(post_process_bind_pose_identity) {
    auto spec = bromesh::builtinHumanoidSpec();
    auto lm = makeHumanoidLandmarks();
    auto mesh = makeSyntheticHumanoid();
    bromesh::VoxelBindOptions opts; opts.maxResolution = 48; opts.smoothIterations = 6;
    auto r = bromesh::autoRig(mesh, spec, lm, opts);

    auto pose = bromesh::bindPose(r.skeleton);
    std::vector<float> world;
    bromesh::computeWorldMatrices(r.skeleton, pose, world);
    auto skinned = mesh;
    bromesh::applySkinning(skinned, r.skin, world.data());
    float maxDelta = 0.0f;
    for (size_t i = 0; i < mesh.positions.size(); ++i) {
        float d = std::fabs(skinned.positions[i] - mesh.positions[i]);
        if (d > maxDelta) maxDelta = d;
    }
    ASSERT(maxDelta < 1e-3f, "bind pose is identity even after smoothing");
}

TEST(post_process_deterministic) {
    auto spec = bromesh::builtinHumanoidSpec();
    auto lm = makeHumanoidLandmarks();
    auto mesh = makeSyntheticHumanoid();
    bromesh::VoxelBindOptions opts; opts.maxResolution = 48; opts.smoothIterations = 4;
    auto a = bromesh::autoRig(mesh, spec, lm, opts);
    auto b = bromesh::autoRig(mesh, spec, lm, opts);
    ASSERT(a.skin.boneWeights == b.skin.boneWeights, "deterministic weights");
    ASSERT(a.skin.boneIndices == b.skin.boneIndices, "deterministic indices");
}

TEST(post_process_side_affinity_intact) {
    // Smoothing must not destroy left/right separation — if it did, posing
    // one arm would drag the other. Re-use the side-affinity check with
    // smoothing on.
    auto spec = bromesh::builtinHumanoidSpec();
    auto lm = makeHumanoidLandmarks();
    auto mesh = makeSyntheticHumanoid();
    bromesh::VoxelBindOptions opts; opts.maxResolution = 48; opts.smoothIterations = 6;
    auto r = bromesh::autoRig(mesh, spec, lm, opts);
    int leftArmBones[3] = { -1, -1, -1 };
    for (size_t i = 0; i < r.skeleton.bones.size(); ++i) {
        const auto& n = r.skeleton.bones[i].name;
        if (n == "upper_arm_L")  leftArmBones[0] = (int)i;
        if (n == "forearm_L")    leftArmBones[1] = (int)i;
        if (n == "shoulder_L")   leftArmBones[2] = (int)i;
    }
    int inspected = 0, withLeft = 0;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float x = mesh.positions[v*3];
        float y = mesh.positions[v*3+1];
        if (x > -0.4f) continue;
        if (y < 0.35f || y > 0.55f) continue;
        ++inspected;
        for (int k = 0; k < 4; ++k) {
            uint32_t bi = r.skin.boneIndices[v*4+k];
            float bw = r.skin.boneWeights[v*4+k];
            if (bw <= 0.0f) continue;
            if ((int)bi == leftArmBones[0] || (int)bi == leftArmBones[1] ||
                (int)bi == leftArmBones[2]) { ++withLeft; break; }
        }
    }
    ASSERT(inspected > 0, "have left-arm verts");
    ASSERT(withLeft * 2 >= inspected, "left-arm verts still bind left-side after smoothing");
}

// --- Phase-5: heuristic landmark detection --------------------------------

TEST(detect_landmarks_humanoid_completeness) {
    auto mesh = makeSyntheticHumanoid();
    auto lm = bromesh::detectHumanoidLandmarks(mesh);
    auto spec = bromesh::builtinHumanoidSpec();
    auto missing = bromesh::missingLandmarks(spec, lm);
    ASSERT(missing.empty(), "all 18 humanoid landmarks detected");
}

TEST(detect_landmarks_humanoid_near_reference) {
    auto mesh = makeSyntheticHumanoid();
    auto detected = bromesh::detectHumanoidLandmarks(mesh);
    auto reference = makeHumanoidLandmarks();

    // Body height of the synthetic mesh is ~1.72. Accept landmarks within
    // 15% of body height from the hand-authored reference — the detector is
    // a geometric heuristic, not a pixel-perfect match.
    auto bbox = bromesh::computeBBox(mesh);
    float H = bbox.extentY();
    float tol = 0.15f * H;

    int checked = 0;
    for (const auto& [name, ref] : reference.points) {
        if (!detected.has(name)) continue;
        auto d = detected.points[name];
        float dx = d[0]-ref[0], dy = d[1]-ref[1], dz = d[2]-ref[2];
        float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
        ASSERT(dist < tol, name.c_str());
        ++checked;
    }
    ASSERT(checked == 18, "checked all 18 landmarks");
}

TEST(detect_landmarks_humanoid_end_to_end) {
    auto mesh = makeSyntheticHumanoid();
    auto spec = bromesh::builtinHumanoidSpec();
    auto lm = bromesh::detectHumanoidLandmarks(mesh);

    bromesh::VoxelBindOptions opts; opts.maxResolution = 48;
    auto r = bromesh::autoRig(mesh, spec, lm, opts);
    ASSERT(r.missingLandmarks.empty(), "no missing landmarks");
    ASSERT(r.skeleton.bones.size() == spec.bones.size(), "bone count");

    // Bind-pose skinning is identity.
    auto pose = bromesh::bindPose(r.skeleton);
    std::vector<float> world;
    bromesh::computeWorldMatrices(r.skeleton, pose, world);
    auto skinned = mesh;
    bromesh::applySkinning(skinned, r.skin, world.data());
    float maxDelta = 0.0f;
    for (size_t i = 0; i < mesh.positions.size(); ++i) {
        float d = std::fabs(skinned.positions[i] - mesh.positions[i]);
        if (d > maxDelta) maxDelta = d;
    }
    ASSERT(maxDelta < 1e-3f, "bind pose skinning is identity");
}

TEST(detect_landmarks_humanoid_deterministic) {
    auto mesh = makeSyntheticHumanoid();
    auto a = bromesh::detectHumanoidLandmarks(mesh);
    auto b = bromesh::detectHumanoidLandmarks(mesh);
    ASSERT(a.points.size() == b.points.size(), "same landmark count");
    for (const auto& [name, pa] : a.points) {
        ASSERT(b.has(name), name.c_str());
        auto pb = b.points[name];
        ASSERT(pa[0] == pb[0] && pa[1] == pb[1] && pa[2] == pb[2],
               "same position both calls");
    }
}

TEST(auto_rig_side_affinity) {
    // Vertices on the left side of the torso should have a left-side bone
    // (upper_arm_L, forearm_L, shoulder_L, or similar) in their top 4.
    auto spec = bromesh::builtinHumanoidSpec();
    auto lm = makeHumanoidLandmarks();
    auto mesh = makeSyntheticHumanoid();
    bromesh::VoxelBindOptions opts; opts.maxResolution = 48;
    auto r = bromesh::autoRig(mesh, spec, lm, opts);

    // Find the left forearm bone index.
    int leftArmBones[3] = { -1, -1, -1 };
    for (size_t i = 0; i < r.skeleton.bones.size(); ++i) {
        const auto& n = r.skeleton.bones[i].name;
        if (n == "upper_arm_L")  leftArmBones[0] = (int)i;
        if (n == "forearm_L")    leftArmBones[1] = (int)i;
        if (n == "shoulder_L")   leftArmBones[2] = (int)i;
    }
    ASSERT(leftArmBones[0] >= 0 && leftArmBones[1] >= 0, "left arm bones present");

    // Count vertices at the extreme -X (left wrist region) that have any
    // left-arm bone in their top 4.
    int inspected = 0, withLeft = 0;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float x = mesh.positions[v * 3];
        float y = mesh.positions[v * 3 + 1];
        if (x > -0.4f) continue;           // only far-left vertices
        if (y < 0.35f || y > 0.55f) continue; // arm band
        ++inspected;
        for (int k = 0; k < 4; ++k) {
            uint32_t bi = r.skin.boneIndices[v * 4 + k];
            float    bw = r.skin.boneWeights[v * 4 + k];
            if (bw <= 0.0f) continue;
            if ((int)bi == leftArmBones[0] ||
                (int)bi == leftArmBones[1] ||
                (int)bi == leftArmBones[2]) { ++withLeft; break; }
        }
    }
    ASSERT(inspected > 0, "inspected some left-arm vertices");
    // Require majority to bind to a left-side bone. Not 100% because box
    // corners touch the torso region.
    ASSERT(withLeft * 2 >= inspected, "majority of left-arm vertices bind to left-side bones");
}
