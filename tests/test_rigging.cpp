#include "test_framework.h"

#include <algorithm>
#include <cmath>

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
#if BROMESH_HAS_GLTF
    ASSERT(!js.empty(), "rig spec JSON should not be empty");
    auto parsed = bromesh::parseRigSpecJSON(js);
    ASSERT(parsed.name == spec.name, "json roundtrip name");
    ASSERT(parsed.bones.size() == spec.bones.size(), "json roundtrip bone count");
    ASSERT(parsed.landmarks.size() == spec.landmarks.size(), "json roundtrip landmark count");
    ASSERT(parsed.sockets.size() == spec.sockets.size(), "json roundtrip socket count");
#else
    ASSERT(js.empty(), "rig spec JSON empty when tinygltf not available");
#endif
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

TEST(skeleton_fit_humanoid_bone_positions) {
    auto spec = bromesh::builtinHumanoidSpec();
    auto lm = makeHumanoidLandmarks();
    auto mesh = makeSyntheticHumanoid();
    auto skel = bromesh::fitSkeleton(spec, lm, mesh);

    auto pose = bromesh::bindPose(skel);
    std::vector<float> world;
    bromesh::computeWorldMatrices(skel, pose, world);

    auto getBonePos = [&](const std::string& name) -> std::array<float, 3> {
        int idx = skel.findBone(name);
        if (idx < 0) return {9999.0f, 9999.0f, 9999.0f};
        return { world[idx * 16 + 12], world[idx * 16 + 13], world[idx * 16 + 14] };
    };

    auto assertClose = [](std::array<float, 3> actual, std::array<float, 3> expected, float tol, const char* msg) {
        float dx = std::fabs(actual[0] - expected[0]);
        float dy = std::fabs(actual[1] - expected[1]);
        float dz = std::fabs(actual[2] - expected[2]);
        ASSERT(dx <= tol && dy <= tol && dz <= tol, msg);
    };

    // wrist_L -> hand_L bone head is landmark:wrist_L (-0.50, 0.45, 0.0)
    assertClose(getBonePos("hand_L"), {-0.50f, 0.45f, 0.00f}, 1e-4f, "hand_L bone at wrist_L landmark");
    // wrist_R -> hand_R bone head is landmark:wrist_R (0.50, 0.45, 0.0)
    assertClose(getBonePos("hand_R"), {0.50f, 0.45f, 0.00f}, 1e-4f, "hand_R bone at wrist_R landmark");
    // ankle_L -> foot_L bone head is landmark:ankle_L (-0.09, -0.88, 0.0)
    assertClose(getBonePos("foot_L"), {-0.09f, -0.88f, 0.00f}, 1e-4f, "foot_L bone at ankle_L landmark");
    // ankle_R -> foot_R bone head is landmark:ankle_R (0.09, -0.88, 0.0)
    assertClose(getBonePos("foot_R"), {0.09f, -0.88f, 0.00f}, 1e-4f, "foot_R bone at ankle_R landmark");
    // head bone head is lerp:neck_base,crown,0.5 -> (0.0, 0.5*(0.55+0.80), 0.0) = (0.0, 0.675, 0.0)
    assertClose(getBonePos("head"), {0.00f, 0.675f, 0.00f}, 1e-4f, "head bone at lerp(neck_base, crown, 0.5)");
    // upper_arm_L head is landmark:shoulder_L (-0.18, 0.50, 0.0)
    assertClose(getBonePos("upper_arm_L"), {-0.18f, 0.50f, 0.00f}, 1e-4f, "upper_arm_L bone at shoulder_L");
    // forearm_L head is landmark:elbow_L (-0.35, 0.45, 0.0)
    assertClose(getBonePos("forearm_L"), {-0.35f, 0.45f, 0.00f}, 1e-4f, "forearm_L bone at elbow_L");
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

    // Bind-pose skinning must leave positions unchanged.
    auto pose = bromesh::bindPose(r.skeleton);
    std::vector<float> joints;
    bromesh::computeSkinningMatrices(r.skeleton, pose, joints);
    auto beforeMesh = mesh;
    bromesh::applySkinning(beforeMesh, r.skin, joints.data());
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

// --- Phase 4 remainder: bone-heat + BBW + weighting dispatch ---------------

// Build a simple manifold capsule-ish mesh by subdividing a box. Used as
// input to the manifold-only weighting paths.
static bromesh::MeshData makeManifoldCapsule() {
    auto m = bromesh::box(0.1f, 0.5f, 0.1f);
    m = bromesh::subdivideMidpoint(m, 2); // smoother discretization
    return m;
}

// Build a simple 2-bone skeleton along Y so per-vertex weights will vary.
static bromesh::Skeleton makeTwoBoneSkeleton() {
    bromesh::Skeleton s;
    bromesh::Bone root; root.name = "root"; root.parent = -1;
    // inverseBind = translate(-0, -0.5, 0) to place the bone at world y=+0.5? We
    // want head in world space at y=-0.5 and tip at y=+0.5. inverseBind is
    // the inverse of world: if world is translate(0, -0.5, 0), inverseBind
    // is translate(0, +0.5, 0). Column-major mat4.
    float ib0[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,+0.5f,0,1};
    std::memcpy(root.inverseBind, ib0, sizeof(ib0));
    bromesh::Bone tip; tip.name = "tip"; tip.parent = 0;
    tip.localT[0] = 0; tip.localT[1] = 1.0f; tip.localT[2] = 0;
    float ib1[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,-0.5f,0,1};
    std::memcpy(tip.inverseBind, ib1, sizeof(ib1));
    s.bones.push_back(root);
    s.bones.push_back(tip);
    return s;
}

TEST(mesh_laplacian_row_sum_zero) {
    // (Lf)[i] = sum over neighbors of w_ij (f[j] - f[i]) — so for f constant,
    // (Lf) must be zero. Equivalent: each row of L sums to zero.
    auto m = makeManifoldCapsule();
    bromesh::SparseCsr L;
    std::vector<double> mass;
    bromesh::assembleCotangentLaplacian(m, L, mass);
    double maxRowSum = 0.0;
    for (int i = 0; i < L.rows; ++i) {
        double s = 0.0;
        for (int k = L.rowStart[i]; k < L.rowStart[i+1]; ++k) s += L.values[k];
        if (std::fabs(s) > maxRowSum) maxRowSum = std::fabs(s);
    }
    ASSERT(maxRowSum < 1e-9, "cot-Laplacian row sums ~ 0");
    double totalMass = 0.0;
    for (double x : mass) totalMass += x;
    ASSERT(totalMass > 1e-6, "mass matrix is positive");
}

TEST(bone_heat_weights_valid) {
    auto mesh = makeManifoldCapsule();
    auto skel = makeTwoBoneSkeleton();
    bromesh::BoneHeatOptions opts;
    auto skin = bromesh::boneHeatWeights(mesh, skel, opts);

    ASSERT(skin.boneCount == 2, "boneCount");
    ASSERT(skin.boneWeights.size() == mesh.vertexCount() * 4, "weights sized");
    // Every vertex: weights sum to ~1 and are non-negative.
    size_t bad = 0;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float s = 0.0f;
        for (int k = 0; k < 4; ++k) {
            float w = skin.boneWeights[v * 4 + k];
            if (!(w >= 0.0f) || w != w) { ++bad; break; }
            s += w;
        }
        if (std::fabs(s - 1.0f) > 1e-3f) ++bad;
    }
    ASSERT(bad == 0, "bone-heat weights normalized and non-negative");
}

TEST(bone_heat_deterministic) {
    auto mesh = makeManifoldCapsule();
    auto skel = makeTwoBoneSkeleton();
    auto a = bromesh::boneHeatWeights(mesh, skel);
    auto b = bromesh::boneHeatWeights(mesh, skel);
    ASSERT(a.boneWeights == b.boneWeights, "bone-heat deterministic");
    ASSERT(a.boneIndices == b.boneIndices, "bone-heat indices deterministic");
}

TEST(bone_heat_spatial_gradient) {
    auto mesh = bromesh::box(0.1f, 1.0f, 0.1f);
    bromesh::translateMesh(mesh, 0.0f, 1.0f, 0.0f); // y in [0, 2]
    mesh = bromesh::weldVertices(mesh, 1e-4f);
    mesh = bromesh::subdivideMidpoint(mesh, 3); // vertices at y = 0, 0.25, 0.5, 0.75, 1.0, ...

    bromesh::Skeleton skel;
    bromesh::Bone b0; b0.name = "bone0"; b0.parent = -1;
    // inverseBind is identity -> head at (0, 0, 0)
    float ib0[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    std::memcpy(b0.inverseBind, ib0, sizeof(ib0));

    bromesh::Bone b1; b1.name = "bone1"; b1.parent = 0;
    b1.localT[0] = 0; b1.localT[1] = 1.0f; b1.localT[2] = 0;
    // inverseBind translate(0, -1, 0) -> head at (0, 1, 0)
    float ib1[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,-1.0f,0,1};
    std::memcpy(b1.inverseBind, ib1, sizeof(ib1));

    skel.bones.push_back(b0);
    skel.bones.push_back(b1);

    bromesh::BoneHeatOptions opts;
    opts.heatStrength = 20.0f;
    auto skin = bromesh::boneHeatWeights(mesh, skel, opts);
    ASSERT(skin.boneCount == 2, "bone_heat_gradient: 2 bones");

    // Check vertices near y = 0.2 (proximal to bone 0, far from bone 1)
    // and near y = 1.8 (proximal to bone 1, far from bone 0)
    int checkedNear0 = 0;
    int checkedNear1 = 0;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float y = mesh.positions[v * 3 + 1];
        float w0 = 0.0f, w1 = 0.0f;
        for (int k = 0; k < 4; ++k) {
            uint32_t bi = skin.boneIndices[v * 4 + k];
            float w = skin.boneWeights[v * 4 + k];
            if (bi == 0) w0 += w;
            if (bi == 1) w1 += w;
        }
        if (y <= 0.25f) {
            checkedNear0++;
            ASSERT(w0 > 0.85f, "bone_heat_gradient: y near 0.2 has bone 0 weight > 0.85");
            ASSERT(w1 < 0.15f, "bone_heat_gradient: y near 0.2 has bone 1 weight < 0.15");
        }
        if (y >= 1.75f) {
            checkedNear1++;
            ASSERT(w1 > 0.85f, "bone_heat_gradient: y near 1.8 has bone 1 weight > 0.85");
            ASSERT(w0 < 0.15f, "bone_heat_gradient: y near 1.8 has bone 0 weight < 0.15");
        }
    }
    ASSERT(checkedNear0 > 0, "bone_heat_gradient: checked vertices near y=0.2");
    ASSERT(checkedNear1 > 0, "bone_heat_gradient: checked vertices near y=1.8");
}

TEST(bbw_weights_valid) {
    auto mesh = makeManifoldCapsule();
    auto skel = makeTwoBoneSkeleton();
    bromesh::BBWOptions opts;
    opts.anchorsPerBone = 2;
    opts.maxIter = 3000;
#if BROMESH_HAS_OSQP
    auto skin = bromesh::bbwWeights(mesh, skel, opts);

    ASSERT(skin.boneCount == 2, "bbw boneCount");
    ASSERT(skin.boneWeights.size() == mesh.vertexCount() * 4, "bbw weights sized");
    size_t bad = 0;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        float s = 0.0f;
        for (int k = 0; k < 4; ++k) {
            float w = skin.boneWeights[v * 4 + k];
            if (!(w >= 0.0f) || w != w) { ++bad; break; }
            s += w;
        }
        if (std::fabs(s - 1.0f) > 1e-2f) ++bad;
    }
    ASSERT(bad == 0, "BBW weights normalized and non-negative");
#else
    auto skin = bromesh::bbwWeights(mesh, skel, opts);
    ASSERT(skin.boneCount == 0, "BBW without OSQP returns empty boneCount");
    ASSERT(skin.boneWeights.empty(), "BBW without OSQP returns empty boneWeights");
#endif
}

TEST(weighting_auto_select_manifold) {
    auto mesh = makeManifoldCapsule();
    auto sel = bromesh::autoSelectWeightingMethod(mesh);
    ASSERT(sel == bromesh::WeightingMethod::BoneHeat, "manifold → bone heat");
}

TEST(weighting_auto_select_non_manifold) {
    // Take a manifold box and punch a hole (remove one triangle) — gives a
    // boundary edge shared by only 1 triangle.
    auto mesh = bromesh::box(0.1f, 0.5f, 0.1f);
    ASSERT(!mesh.indices.empty(), "sanity: box has triangles");
    mesh.indices.resize(mesh.indices.size() - 3);
    auto sel = bromesh::autoSelectWeightingMethod(mesh);
    ASSERT(sel == bromesh::WeightingMethod::VoxelBind, "non-manifold → voxel");
}

TEST(builtin_rig_spec_dispatcher) {
    ASSERT(!bromesh::builtinRigSpec("humanoid").bones.empty(),  "humanoid spec");
    ASSERT(!bromesh::builtinRigSpec("quadruped").bones.empty(), "quadruped spec");
    ASSERT(!bromesh::builtinRigSpec("hexapod").bones.empty(),   "hexapod spec");
    ASSERT(!bromesh::builtinRigSpec("octopod").bones.empty(),   "octopod spec");
    ASSERT(bromesh::builtinRigSpec("nonsense").bones.empty(),   "unknown → empty");
}

TEST(weighting_method_name_roundtrip) {
    using bromesh::WeightingMethod;
    for (auto m : { WeightingMethod::Auto, WeightingMethod::VoxelBind,
                    WeightingMethod::BoneHeat, WeightingMethod::BBW }) {
        ASSERT(bromesh::parseWeightingMethod(bromesh::weightingMethodName(m)) == m,
               "method name round-trips");
    }
    ASSERT(bromesh::parseWeightingMethod("garbage") == WeightingMethod::Auto,
           "unknown method → auto");
}

TEST(validate_skin_matches_test_bar) {
    auto spec = bromesh::builtinHumanoidSpec();
    auto lm = makeHumanoidLandmarks();
    auto mesh = makeSyntheticHumanoid();
    bromesh::VoxelBindOptions opts; opts.maxResolution = 48;
    auto r = bromesh::autoRig(mesh, spec, lm, opts);
    auto v = bromesh::validateSkin(mesh, r.skin, 4);
    ASSERT(v.vertexCount == mesh.vertexCount(),      "vertexCount");
    ASSERT(v.clean(),                                 "validateSkin clean on good rig");
    ASSERT(v.maxInfluencesObserved <= 4,             "at most 4 influences");
}

TEST(auto_rig_new_options_path_compat) {
    // Back-compat overload (VoxelBindOptions) still works and matches the
    // new overload when method is forced to VoxelBind.
    auto spec = bromesh::builtinHumanoidSpec();
    auto lm = makeHumanoidLandmarks();
    auto mesh = makeSyntheticHumanoid();
    bromesh::VoxelBindOptions bOpts; bOpts.maxResolution = 48;
    auto legacy = bromesh::autoRig(mesh, spec, lm, bOpts);

    bromesh::WeightingOptions wOpts;
    wOpts.method = bromesh::WeightingMethod::VoxelBind;
    wOpts.voxel = bOpts;
    wOpts.smoothIterations = bOpts.smoothIterations;
    wOpts.smoothAlpha = bOpts.smoothAlpha;
    wOpts.minWeight = bOpts.minWeight;
    auto newPath = bromesh::autoRig(mesh, spec, lm, wOpts);

    ASSERT(legacy.skin.boneWeights == newPath.skin.boneWeights,
           "new/legacy autoRig produce same weights");
    ASSERT(legacy.methodUsed == bromesh::WeightingMethod::VoxelBind,
           "legacy path reports VoxelBind");
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
    std::vector<float> joints;
    bromesh::computeSkinningMatrices(r.skeleton, pose, joints);
    auto skinned = mesh;
    bromesh::applySkinning(skinned, r.skin, joints.data());
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

    auto bbox = bromesh::computeBBox(mesh);
    float H = bromath::aextent(bbox).y;
    float tol = 0.05f * H;

    int checked = 0;
    for (const auto& [name, ref] : reference.points) {
        if (!detected.has(name)) continue;
        auto d = detected.points.at(name);
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
    std::vector<float> joints;
    bromesh::computeSkinningMatrices(r.skeleton, pose, joints);
    auto skinned = mesh;
    bromesh::applySkinning(skinned, r.skin, joints.data());
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
