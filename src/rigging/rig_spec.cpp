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

// ---- builtinQuadrupedSpec --------------------------------------------------

RigSpec builtinQuadrupedSpec() {
    RigSpec s;
    s.name = "quadruped";
    s.symmetric = true;

    auto lm = [&](const char* name, const char* mirror = "") {
        RigSpec::LandmarkDecl d;
        d.name = name;
        d.mirror = mirror;
        s.landmarks.push_back(d);
    };
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

    // Axial landmarks. +Z is forward (snout direction); tail extends -Z.
    lm("pelvis");
    lm("chest");
    lm("neck_base");
    lm("crown");
    lm("muzzle");
    lm("tail_base");
    lm("tail_tip");
    // Front legs.
    lm("fshoulder_L", "fshoulder_R");
    lm("fshoulder_R", "fshoulder_L");
    lm("felbow_L",    "felbow_R");
    lm("felbow_R",    "felbow_L");
    lm("fpaw_L",      "fpaw_R");
    lm("fpaw_R",      "fpaw_L");
    // Hind legs.
    lm("hip_L",       "hip_R");
    lm("hip_R",       "hip_L");
    lm("hknee_L",     "hknee_R");
    lm("hknee_R",     "hknee_L");
    lm("hpaw_L",      "hpaw_R");
    lm("hpaw_R",      "hpaw_L");

    // Axial chain: pelvis -> spine_01 -> chest -> neck -> head.
    bone("pelvis",   "",         "landmark:pelvis",                 "lerp:pelvis,chest,0.33");
    bone("spine_01", "pelvis",   "lerp:pelvis,chest,0.33",          "lerp:pelvis,chest,0.66");
    bone("chest",    "spine_01", "lerp:pelvis,chest,0.66",          "landmark:neck_base");
    bone("neck",     "chest",    "landmark:neck_base",              "lerp:neck_base,crown,0.5");
    bone("head",     "neck",     "lerp:neck_base,crown,0.5",        "landmark:muzzle");

    // Tail chain hangs off pelvis.
    bone("tail_01", "pelvis",  "landmark:tail_base",                "lerp:tail_base,tail_tip,0.33");
    bone("tail_02", "tail_01", "lerp:tail_base,tail_tip,0.33",      "lerp:tail_base,tail_tip,0.66");
    bone("tail_03", "tail_02", "lerp:tail_base,tail_tip,0.66",      "landmark:tail_tip");

    // Front legs (parent chest).
    bone("fl_upper_L",  "chest",      "landmark:fshoulder_L",  "landmark:felbow_L");
    bone("fl_upper_R",  "chest",      "landmark:fshoulder_R",  "landmark:felbow_R");
    bone("fl_lower_L",  "fl_upper_L", "landmark:felbow_L",     "landmark:fpaw_L");
    bone("fl_lower_R",  "fl_upper_R", "landmark:felbow_R",     "landmark:fpaw_R");
    bone("fl_paw_L",    "fl_lower_L", "landmark:fpaw_L",       "offset:fpaw_L,0,0,0.06", false);
    bone("fl_paw_R",    "fl_lower_R", "landmark:fpaw_R",       "offset:fpaw_R,0,0,0.06", false);

    // Hind legs (parent pelvis).
    bone("hl_upper_L",  "pelvis",     "landmark:hip_L",        "landmark:hknee_L");
    bone("hl_upper_R",  "pelvis",     "landmark:hip_R",        "landmark:hknee_R");
    bone("hl_lower_L",  "hl_upper_L", "landmark:hknee_L",      "landmark:hpaw_L");
    bone("hl_lower_R",  "hl_upper_R", "landmark:hknee_R",      "landmark:hpaw_R");
    bone("hl_paw_L",    "hl_lower_L", "landmark:hpaw_L",       "offset:hpaw_L,0,0,0.06", false);
    bone("hl_paw_R",    "hl_lower_R", "landmark:hpaw_R",       "offset:hpaw_R,0,0,0.06", false);

    // IK chain tips (paws).
    for (auto& b : s.bones) {
        if (b.name == "fl_paw_L" || b.name == "fl_paw_R" ||
            b.name == "hl_paw_L" || b.name == "hl_paw_R") {
            b.ikChain = true;
        }
    }

    auto socket = [&](const char* name, const char* boneName) {
        RigSpec::SocketDecl sk;
        sk.name = name;
        sk.bone = boneName;
        s.sockets.push_back(sk);
    };
    socket("head_attach", "head");
    socket("saddle",      "chest");

    return s;
}

// ---- builtinHexapodSpec ----------------------------------------------------

RigSpec builtinHexapodSpec() {
    RigSpec s;
    s.name = "hexapod";
    s.symmetric = true;

    auto lm = [&](const char* name, const char* mirror = "") {
        RigSpec::LandmarkDecl d;
        d.name = name;
        d.mirror = mirror;
        s.landmarks.push_back(d);
    };
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

    // Body landmarks. +Z forward (head direction); abdomen trails at -Z.
    lm("abdomen_tip");
    lm("abdomen");
    lm("thorax");
    lm("head");
    lm("head_tip");
    // Three bilateral leg pairs. Three segments per leg (coxa / femur / tibia
    // tip), represented by three landmarks per leg side.
    const char* prefixes[3] = { "front", "mid", "rear" };
    for (const char* p : prefixes) {
        std::string hipL = std::string(p) + "_hip_L";
        std::string hipR = std::string(p) + "_hip_R";
        std::string kneeL = std::string(p) + "_knee_L";
        std::string kneeR = std::string(p) + "_knee_R";
        std::string footL = std::string(p) + "_foot_L";
        std::string footR = std::string(p) + "_foot_R";
        lm(hipL.c_str(),  hipR.c_str());
        lm(hipR.c_str(),  hipL.c_str());
        lm(kneeL.c_str(), kneeR.c_str());
        lm(kneeR.c_str(), kneeL.c_str());
        lm(footL.c_str(), footR.c_str());
        lm(footR.c_str(), footL.c_str());
    }

    // Body chain: abdomen -> thorax -> head.
    bone("abdomen", "",        "landmark:abdomen",           "landmark:thorax");
    bone("thorax",  "abdomen", "landmark:thorax",            "landmark:head");
    bone("head",    "thorax",  "landmark:head",              "landmark:head_tip");

    // Legs all parent on thorax.
    for (const char* p : prefixes) {
        std::string coxaL  = std::string(p) + "_coxa_L";
        std::string coxaR  = std::string(p) + "_coxa_R";
        std::string femurL = std::string(p) + "_femur_L";
        std::string femurR = std::string(p) + "_femur_R";
        std::string tibiaL = std::string(p) + "_tibia_L";
        std::string tibiaR = std::string(p) + "_tibia_R";
        std::string hipL   = std::string("landmark:") + p + "_hip_L";
        std::string hipR   = std::string("landmark:") + p + "_hip_R";
        std::string kneeL  = std::string("landmark:") + p + "_knee_L";
        std::string kneeR  = std::string("landmark:") + p + "_knee_R";
        std::string footL  = std::string("landmark:") + p + "_foot_L";
        std::string footR  = std::string("landmark:") + p + "_foot_R";
        std::string footLtail = std::string("offset:") + p + "_foot_L,0,-0.02,0";
        std::string footRtail = std::string("offset:") + p + "_foot_R,0,-0.02,0";

        bone(coxaL.c_str(),  "thorax",         hipL.c_str(),  kneeL.c_str());
        bone(coxaR.c_str(),  "thorax",         hipR.c_str(),  kneeR.c_str());
        bone(femurL.c_str(), coxaL.c_str(),    kneeL.c_str(), footL.c_str());
        bone(femurR.c_str(), coxaR.c_str(),    kneeR.c_str(), footR.c_str());
        bone(tibiaL.c_str(), femurL.c_str(),   footL.c_str(), footLtail.c_str(), false);
        bone(tibiaR.c_str(), femurR.c_str(),   footR.c_str(), footRtail.c_str(), false);
    }

    // IK chain tips (tibias — final leg segment).
    for (auto& b : s.bones) {
        if (b.name.size() > 7 && b.name.find("_tibia_") != std::string::npos) {
            b.ikChain = true;
        }
    }

    RigSpec::SocketDecl head_attach;
    head_attach.name = "head_attach";
    head_attach.bone = "head";
    s.sockets.push_back(head_attach);

    return s;
}

// ---- builtinOctopodSpec ----------------------------------------------------

RigSpec builtinOctopodSpec() {
    RigSpec s;
    s.name = "octopod";
    s.symmetric = true;

    auto lm = [&](const char* name, const char* mirror = "") {
        RigSpec::LandmarkDecl d;
        d.name = name;
        d.mirror = mirror;
        s.landmarks.push_back(d);
    };
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

    lm("body");
    lm("head");
    // Four bilateral arm pairs. Arms labelled a..d from front to back.
    // Two segments per arm; three landmarks per arm side.
    const char* arms[4] = { "arm_a", "arm_b", "arm_c", "arm_d" };
    for (const char* a : arms) {
        std::string hipL  = std::string(a) + "_hip_L";
        std::string hipR  = std::string(a) + "_hip_R";
        std::string midL  = std::string(a) + "_mid_L";
        std::string midR  = std::string(a) + "_mid_R";
        std::string tipL  = std::string(a) + "_tip_L";
        std::string tipR  = std::string(a) + "_tip_R";
        lm(hipL.c_str(), hipR.c_str());
        lm(hipR.c_str(), hipL.c_str());
        lm(midL.c_str(), midR.c_str());
        lm(midR.c_str(), midL.c_str());
        lm(tipL.c_str(), tipR.c_str());
        lm(tipR.c_str(), tipL.c_str());
    }

    // Body chain.
    bone("body", "",     "landmark:body", "landmark:head");
    bone("head", "body", "landmark:head", "offset:head,0,0.05,0");

    // Arms all parent on body.
    for (const char* a : arms) {
        std::string upperL = std::string(a) + "_upper_L";
        std::string upperR = std::string(a) + "_upper_R";
        std::string lowerL = std::string(a) + "_lower_L";
        std::string lowerR = std::string(a) + "_lower_R";
        std::string hipL   = std::string("landmark:") + a + "_hip_L";
        std::string hipR   = std::string("landmark:") + a + "_hip_R";
        std::string midL   = std::string("landmark:") + a + "_mid_L";
        std::string midR   = std::string("landmark:") + a + "_mid_R";
        std::string tipL   = std::string("landmark:") + a + "_tip_L";
        std::string tipR   = std::string("landmark:") + a + "_tip_R";

        bone(upperL.c_str(), "body",          hipL.c_str(), midL.c_str());
        bone(upperR.c_str(), "body",          hipR.c_str(), midR.c_str());
        bone(lowerL.c_str(), upperL.c_str(),  midL.c_str(), tipL.c_str(), false);
        bone(lowerR.c_str(), upperR.c_str(),  midR.c_str(), tipR.c_str(), false);
    }

    // IK chain tips (lower arm segments).
    for (auto& b : s.bones) {
        if (b.name.size() > 2 && b.name.find("_lower_") != std::string::npos) {
            b.ikChain = true;
        }
    }

    RigSpec::SocketDecl body_attach;
    body_attach.name = "body_attach";
    body_attach.bone = "body";
    s.sockets.push_back(body_attach);

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

RigSpec builtinRigSpec(const std::string& name) {
    if (name == "humanoid")  return builtinHumanoidSpec();
    if (name == "quadruped") return builtinQuadrupedSpec();
    if (name == "hexapod")   return builtinHexapodSpec();
    if (name == "octopod")   return builtinOctopodSpec();
    return {};
}

} // namespace bromesh
