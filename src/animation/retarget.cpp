#include "bromesh/animation/retarget.h"

#include <cstring>

namespace bromesh {

Animation retargetAnimation(const Animation& anim,
                             const Skeleton& srcSkeleton,
                             const Skeleton& dstSkeleton) {
    Animation out;
    out.name = anim.name;
    out.duration = anim.duration;

    for (const auto& ch : anim.channels) {
        if (ch.boneIndex < 0 || ch.boneIndex >= (int)srcSkeleton.bones.size()) continue;
        const std::string& name = srcSkeleton.bones[ch.boneIndex].name;
        int dstIdx = dstSkeleton.findBone(name);
        if (dstIdx < 0) continue;
        AnimChannel nc = ch;
        nc.boneIndex = dstIdx;
        out.channels.push_back(std::move(nc));
    }
    return out;
}

int findBoneBySuffix(const Skeleton& skeleton, const std::string& suffix) {
    for (size_t i = 0; i < skeleton.bones.size(); ++i) {
        const std::string& n = skeleton.bones[i].name;
        if (n.size() >= suffix.size() &&
            n.compare(n.size() - suffix.size(), suffix.size(), suffix) == 0) {
            return (int)i;
        }
    }
    return -1;
}

// Look for a bone whose name matches any of the given candidates exactly or
// as a suffix (to handle DEF- / mixamorig: prefixes).
static int findByCandidates(const Skeleton& sk, std::initializer_list<const char*> candidates) {
    for (const char* c : candidates) {
        int i = sk.findBone(c);
        if (i >= 0) return i;
    }
    for (const char* c : candidates) {
        int i = findBoneBySuffix(sk, std::string(c));
        if (i >= 0) return i;
    }
    return -1;
}

int addRigifySockets(Skeleton& skeleton) {
    struct Entry { const char* socketName; std::initializer_list<const char*> names; };
    Entry entries[] = {
        // Rigify deform bone names use ".L"/".R" suffixes; Mixamo uses LeftHand/RightHand.
        { "hand.R",  { "DEF-hand.R",  "hand.R",  "mixamorig:RightHand" } },
        { "hand.L",  { "DEF-hand.L",  "hand.L",  "mixamorig:LeftHand"  } },
        { "foot.R",  { "DEF-foot.R",  "foot.R",  "mixamorig:RightFoot" } },
        { "foot.L",  { "DEF-foot.L",  "foot.L",  "mixamorig:LeftFoot"  } },
        { "head",    { "DEF-head",    "head",    "mixamorig:Head"      } },
        { "chest",   { "DEF-chest",   "DEF-spine.003", "chest", "mixamorig:Spine2" } },
        { "pelvis",  { "DEF-pelvis",  "DEF-spine", "pelvis", "mixamorig:Hips" } },
    };

    int added = 0;
    for (const auto& e : entries) {
        int b = findByCandidates(skeleton, e.names);
        if (b < 0) continue;
        if (skeleton.findSocket(e.socketName) >= 0) continue;
        Socket s;
        s.name = e.socketName;
        s.bone = b;
        // Identity offset
        for (int k = 0; k < 16; ++k) s.offset[k] = 0.0f;
        s.offset[0] = s.offset[5] = s.offset[10] = s.offset[15] = 1.0f;
        skeleton.sockets.push_back(std::move(s));
        ++added;
    }
    return added;
}

} // namespace bromesh
