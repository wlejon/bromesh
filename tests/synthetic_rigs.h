#pragma once

// Shared synthetic mesh + landmark factories used by both the rigging and
// animation test suites. Each creature matches one of the builtin rig specs
// (humanoid, quadruped, hexapod, octopod) and uses the project-wide
// convention: +Y up, +Z forward, +X right.

#include "test_framework.h"

namespace synthrig {

inline bromesh::MeshData boxAt(float cx, float cy, float cz,
                               float hx, float hy, float hz) {
    auto m = bromesh::box(hx, hy, hz);
    bromesh::translateMesh(m, cx, cy, cz);
    return m;
}

// ---- Humanoid ------------------------------------------------------------

inline bromesh::MeshData makeHumanoidMesh() {
    std::vector<bromesh::MeshData> parts;
    parts.push_back(boxAt(0.0f,  0.25f, 0.0f, 0.18f, 0.30f, 0.10f));
    parts.push_back(boxAt(0.0f,  0.70f, 0.0f, 0.12f, 0.12f, 0.12f));
    parts.push_back(boxAt(-0.35f, 0.45f, 0.0f, 0.18f, 0.06f, 0.06f));
    parts.push_back(boxAt( 0.35f, 0.45f, 0.0f, 0.18f, 0.06f, 0.06f));
    parts.push_back(boxAt(-0.09f,-0.45f, 0.0f, 0.06f, 0.45f, 0.06f));
    parts.push_back(boxAt( 0.09f,-0.45f, 0.0f, 0.06f, 0.45f, 0.06f));
    return bromesh::mergeMeshes(parts);
}

inline bromesh::Landmarks makeHumanoidLandmarks() {
    bromesh::Landmarks lm;
    lm.set("pelvis",    0.00f, -0.05f, 0.00f);
    lm.set("chest",     0.00f,  0.45f, 0.00f);
    lm.set("neck_base", 0.00f,  0.55f, 0.00f);
    lm.set("crown",     0.00f,  0.80f, 0.00f);
    lm.set("shoulder_L",-0.18f, 0.50f, 0.00f);
    lm.set("shoulder_R", 0.18f, 0.50f, 0.00f);
    lm.set("elbow_L",   -0.35f, 0.45f, 0.00f);
    lm.set("elbow_R",    0.35f, 0.45f, 0.00f);
    lm.set("wrist_L",   -0.50f, 0.45f, 0.00f);
    lm.set("wrist_R",    0.50f, 0.45f, 0.00f);
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

// ---- Quadruped -----------------------------------------------------------

inline bromesh::MeshData makeQuadrupedMesh() {
    std::vector<bromesh::MeshData> parts;
    parts.push_back(boxAt(0.00f, 0.45f, 0.00f, 0.12f, 0.12f, 0.30f));
    parts.push_back(boxAt(0.00f, 0.55f, 0.40f, 0.08f, 0.08f, 0.10f));
    parts.push_back(boxAt(0.00f, 0.50f,-0.45f, 0.03f, 0.03f, 0.15f));
    parts.push_back(boxAt(-0.08f, 0.18f, 0.20f, 0.04f, 0.27f, 0.04f));
    parts.push_back(boxAt( 0.08f, 0.18f, 0.20f, 0.04f, 0.27f, 0.04f));
    parts.push_back(boxAt(-0.08f, 0.18f,-0.20f, 0.04f, 0.27f, 0.04f));
    parts.push_back(boxAt( 0.08f, 0.18f,-0.20f, 0.04f, 0.27f, 0.04f));
    return bromesh::mergeMeshes(parts);
}

inline bromesh::Landmarks makeQuadrupedLandmarks() {
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

inline bromesh::MeshData makeHexapodMesh() {
    std::vector<bromesh::MeshData> parts;
    parts.push_back(boxAt(0.00f, 0.20f,-0.20f, 0.08f, 0.08f, 0.15f));
    parts.push_back(boxAt(0.00f, 0.20f, 0.00f, 0.10f, 0.10f, 0.12f));
    parts.push_back(boxAt(0.00f, 0.20f, 0.18f, 0.06f, 0.06f, 0.06f));
    for (int side = 0; side < 2; ++side) {
        float sx = (side == 0) ? -1.0f : 1.0f;
        for (int row = 0; row < 3; ++row) {
            float z = 0.08f - row * 0.08f;
            parts.push_back(boxAt(sx * 0.12f, 0.10f, z, 0.03f, 0.15f, 0.03f));
        }
    }
    return bromesh::mergeMeshes(parts);
}

inline bromesh::Landmarks makeHexapodLandmarks() {
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

inline bromesh::MeshData makeOctopodMesh() {
    std::vector<bromesh::MeshData> parts;
    parts.push_back(boxAt(0.00f, 0.25f, 0.00f, 0.14f, 0.08f, 0.14f));
    parts.push_back(boxAt(0.00f, 0.35f, 0.00f, 0.06f, 0.04f, 0.06f));
    float armZ[4] = { 0.15f, 0.05f, -0.05f, -0.15f };
    for (int i = 0; i < 4; ++i) {
        float z = armZ[i];
        parts.push_back(boxAt(-0.22f, 0.15f, z, 0.08f, 0.03f, 0.03f));
        parts.push_back(boxAt( 0.22f, 0.15f, z, 0.08f, 0.03f, 0.03f));
    }
    return bromesh::mergeMeshes(parts);
}

inline bromesh::Landmarks makeOctopodLandmarks() {
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

} // namespace synthrig
