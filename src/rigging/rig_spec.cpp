#include "bromesh/rigging/rig_spec.h"
#include "bromesh/rigging/landmarks.h"

#include <fstream>
#include <sstream>

#if BROMESH_HAS_GLTF
#include "json.hpp"
#endif

namespace bromesh {

// ---- missingLandmarks ------------------------------------------------------

std::vector<std::string> missingLandmarks(const RigSpec& spec,
                                          const Landmarks& lm) {
    std::vector<std::string> missing;
    for (const auto& d : spec.landmarks) {
        if (!lm.has(d.name)) missing.push_back(d.name);
    }
    return missing;
}

// ---- builtinHumanoidSpec ---------------------------------------------------

RigSpec builtinHumanoidSpec() {
    RigSpec s;
    s.name = "humanoid";
    s.symmetric = true;

    auto lm = [&](const char* name, const char* mirror = "") {
        RigSpec::LandmarkDecl d;
        d.name = name;
        d.mirror = mirror;
        s.landmarks.push_back(d);
    };

    // 16 user-marked landmarks. Conventional T-pose; +Y up, +Z forward.
    lm("pelvis");
    lm("chest");
    lm("neck_base");
    lm("crown");
    lm("shoulder_L", "shoulder_R");
    lm("shoulder_R", "shoulder_L");
    lm("elbow_L",    "elbow_R");
    lm("elbow_R",    "elbow_L");
    lm("wrist_L",    "wrist_R");
    lm("wrist_R",    "wrist_L");
    lm("hip_L",      "hip_R");
    lm("hip_R",      "hip_L");
    lm("knee_L",     "knee_R");
    lm("knee_R",     "knee_L");
    lm("ankle_L",    "ankle_R");
    lm("ankle_R",    "ankle_L");
    lm("toe_L",      "toe_R");
    lm("toe_R",      "toe_L");

    auto bone = [&](const char* name, const char* parent,
                    const char* head, const char* tail,
                    bool interior = true) {
        RigSpec::BoneDecl b;
        b.name = name;
        b.parent = parent;
        b.head = head;
        b.tail = tail;
        b.interior = interior;
        s.bones.push_back(b);
    };

    // Spine. 3 spine segments between pelvis and neck_base.
    bone("pelvis",   "",         "landmark:pelvis",                 "landmark:chest");
    bone("spine_01", "pelvis",   "lerp:pelvis,chest,0.33",          "lerp:pelvis,chest,0.66");
    bone("spine_02", "spine_01", "lerp:pelvis,chest,0.66",          "landmark:chest");
    bone("spine_03", "spine_02", "landmark:chest",                  "landmark:neck_base");
    bone("neck",     "spine_03", "landmark:neck_base",              "lerp:neck_base,crown,0.5");
    bone("head",     "neck",     "lerp:neck_base,crown,0.5",        "landmark:crown");

    // Left arm.
    bone("shoulder_L",  "spine_03",   "lerp:neck_base,shoulder_L,0.6", "landmark:shoulder_L");
    bone("upper_arm_L", "shoulder_L", "landmark:shoulder_L",           "landmark:elbow_L");
    bone("forearm_L",   "upper_arm_L","landmark:elbow_L",              "landmark:wrist_L");
    bone("hand_L",      "forearm_L",  "landmark:wrist_L",              "offset:wrist_L,0.08,0,0", false);

    // Right arm.
    bone("shoulder_R",  "spine_03",   "lerp:neck_base,shoulder_R,0.6", "landmark:shoulder_R");
    bone("upper_arm_R", "shoulder_R", "landmark:shoulder_R",           "landmark:elbow_R");
    bone("forearm_R",   "upper_arm_R","landmark:elbow_R",              "landmark:wrist_R");
    bone("hand_R",      "forearm_R",  "landmark:wrist_R",              "offset:wrist_R,-0.08,0,0", false);

    // Left leg.
    bone("upper_leg_L", "pelvis",     "landmark:hip_L",    "landmark:knee_L");
    bone("lower_leg_L", "upper_leg_L","landmark:knee_L",   "landmark:ankle_L");
    bone("foot_L",      "lower_leg_L","landmark:ankle_L",  "landmark:toe_L", false);
    bone("toe_L",       "foot_L",     "landmark:toe_L",    "offset:toe_L,0,0,0.04", false);

    // Right leg.
    bone("upper_leg_R", "pelvis",     "landmark:hip_R",    "landmark:knee_R");
    bone("lower_leg_R", "upper_leg_R","landmark:knee_R",   "landmark:ankle_R");
    bone("foot_R",      "lower_leg_R","landmark:ankle_R",  "landmark:toe_R", false);
    bone("toe_R",       "foot_R",     "landmark:toe_R",    "offset:toe_R,0,0,0.04", false);

    // Mark IK chain tips for future locomotion / IK solvers.
    for (auto& b : s.bones) {
        if (b.name == "hand_L" || b.name == "hand_R" ||
            b.name == "foot_L" || b.name == "foot_R") {
            b.ikChain = true;
        }
    }

    // A couple of default sockets. Identity offset at the specified bone.
    auto socket = [&](const char* name, const char* boneName) {
        RigSpec::SocketDecl sk;
        sk.name = name;
        sk.bone = boneName;
        s.sockets.push_back(sk);
    };
    socket("hand_L_grip", "hand_L");
    socket("hand_R_grip", "hand_R");
    socket("head_attach", "head");

    return s;
}

// ---- JSON parse / serialize ------------------------------------------------

#if BROMESH_HAS_GLTF

namespace {

using json = nlohmann::json;

RigSpec fromJson(const json& j) {
    RigSpec s;
    if (j.contains("name") && j["name"].is_string()) s.name = j["name"].get<std::string>();
    if (j.contains("symmetric") && j["symmetric"].is_boolean()) s.symmetric = j["symmetric"].get<bool>();

    if (j.contains("landmarks") && j["landmarks"].is_array()) {
        for (const auto& e : j["landmarks"]) {
            RigSpec::LandmarkDecl d;
            if (e.contains("name")) d.name = e["name"].get<std::string>();
            if (e.contains("mirror")) d.mirror = e["mirror"].get<std::string>();
            s.landmarks.push_back(std::move(d));
        }
    }
    if (j.contains("bones") && j["bones"].is_array()) {
        for (const auto& e : j["bones"]) {
            RigSpec::BoneDecl b;
            if (e.contains("name"))     b.name = e["name"].get<std::string>();
            if (e.contains("parent"))   b.parent = e["parent"].get<std::string>();
            if (e.contains("head"))     b.head = e["head"].get<std::string>();
            if (e.contains("tail"))     b.tail = e["tail"].get<std::string>();
            if (e.contains("length"))   b.length = e["length"].get<float>();
            if (e.contains("ikChain"))  b.ikChain = e["ikChain"].get<bool>();
            if (e.contains("interior")) b.interior = e["interior"].get<bool>();
            s.bones.push_back(std::move(b));
        }
    }
    if (j.contains("sockets") && j["sockets"].is_array()) {
        for (const auto& e : j["sockets"]) {
            RigSpec::SocketDecl sk;
            if (e.contains("name")) sk.name = e["name"].get<std::string>();
            if (e.contains("bone")) sk.bone = e["bone"].get<std::string>();
            if (e.contains("offset") && e["offset"].is_array() && e["offset"].size() == 16) {
                for (size_t i = 0; i < 16; ++i) sk.offset[i] = e["offset"][i].get<float>();
            }
            s.sockets.push_back(std::move(sk));
        }
    }
    return s;
}

json toJson(const RigSpec& s) {
    json j;
    j["name"] = s.name;
    j["symmetric"] = s.symmetric;
    j["landmarks"] = json::array();
    for (const auto& d : s.landmarks) {
        json e;
        e["name"] = d.name;
        if (!d.mirror.empty()) e["mirror"] = d.mirror;
        j["landmarks"].push_back(std::move(e));
    }
    j["bones"] = json::array();
    for (const auto& b : s.bones) {
        json e;
        e["name"] = b.name;
        e["parent"] = b.parent;
        e["head"] = b.head;
        if (!b.tail.empty()) e["tail"] = b.tail;
        if (b.length != 0.0f) e["length"] = b.length;
        if (b.ikChain) e["ikChain"] = true;
        if (!b.interior) e["interior"] = false;
        j["bones"].push_back(std::move(e));
    }
    j["sockets"] = json::array();
    for (const auto& sk : s.sockets) {
        json e;
        e["name"] = sk.name;
        e["bone"] = sk.bone;
        e["offset"] = json::array();
        for (int i = 0; i < 16; ++i) e["offset"].push_back(sk.offset[i]);
        j["sockets"].push_back(std::move(e));
    }
    return j;
}

} // namespace

RigSpec parseRigSpecJSON(const std::string& jsonText) {
    try {
        json j = json::parse(jsonText);
        return fromJson(j);
    } catch (...) {
        return {};
    }
}

std::string serializeRigSpecJSON(const RigSpec& spec) {
    return toJson(spec).dump(2);
}

#else // !BROMESH_HAS_GLTF

RigSpec parseRigSpecJSON(const std::string&) { return {}; }
std::string serializeRigSpecJSON(const RigSpec&) { return {}; }

#endif

RigSpec loadRigSpecFile(const std::string& path) {
    std::ifstream f(path);
    if (!f) return {};
    std::stringstream buf;
    buf << f.rdbuf();
    return parseRigSpecJSON(buf.str());
}

} // namespace bromesh
