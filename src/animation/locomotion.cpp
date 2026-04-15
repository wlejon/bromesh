#include "bromesh/animation/locomotion.h"
#include "bromesh/animation/pose.h"
#include "bromesh/animation/ik.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>

namespace bromesh {

namespace {

constexpr float kPi = 3.14159265358979323846f;

void allBoneHeadsWorld(const Skeleton& sk, const Pose& pose,
                       std::vector<float>& out) {
    std::vector<float> world;
    computeWorldMatrices(sk, pose, world);
    const size_t n = sk.bones.size();
    out.resize(n * 3);
    for (size_t i = 0; i < n; ++i) {
        out[i*3+0] = world[i*16 + 12];
        out[i*3+1] = world[i*16 + 13];
        out[i*3+2] = world[i*16 + 14];
    }
}

int countChildren(const Skeleton& sk, int parent) {
    int n = 0;
    for (const auto& b : sk.bones)
        if (b.parent == parent) ++n;
    return n;
}

float proj(const float v[3], const float axis[3]) {
    return v[0]*axis[0] + v[1]*axis[1] + v[2]*axis[2];
}

void vnorm(float v[3]) {
    float l = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if (l > 1e-6f) { v[0]/=l; v[1]/=l; v[2]/=l; }
}

} // namespace

// ---- identifyLegChains ----------------------------------------------------

std::vector<LegChain> identifyLegChains(const Skeleton& sk, const RigSpec& spec,
                                        const float upAxisIn[3]) {
    float upAxis[3] = {0.0f, 1.0f, 0.0f};
    if (upAxisIn) { upAxis[0]=upAxisIn[0]; upAxis[1]=upAxisIn[1]; upAxis[2]=upAxisIn[2]; }
    vnorm(upAxis);

    std::vector<LegChain> chains;
    if (sk.bones.empty() || spec.bones.empty()) return chains;

    Pose bind = bindPose(sk);
    std::vector<float> heads;
    allBoneHeadsWorld(sk, bind, heads);

    // Skeleton extent along up axis.
    float minU = std::numeric_limits<float>::max();
    float maxU = -std::numeric_limits<float>::max();
    for (size_t i = 0; i < sk.bones.size(); ++i) {
        float u = proj(&heads[i*3], upAxis);
        minU = std::min(minU, u);
        maxU = std::max(maxU, u);
    }
    const float midU = 0.5f * (minU + maxU);

    for (const auto& decl : spec.bones) {
        if (!decl.ikChain) continue;
        int tipIdx = sk.findBone(decl.name);
        if (tipIdx < 0) continue;

        // Walk up parent links until we hit a fork (multi-child ancestor) or
        // the root. The top of the chain is the fork's direct child (the
        // "hip"); the fork itself is excluded.
        std::vector<int> chain{ tipIdx };
        int cur = tipIdx;
        while (true) {
            int p = sk.bones[cur].parent;
            if (p < 0) break;
            if (countChildren(sk, p) > 1) break;
            chain.push_back(p);
            cur = p;
        }
        std::reverse(chain.begin(), chain.end());

        LegChain lc;
        lc.bones = std::move(chain);
        lc.tipBone = tipIdx;
        lc.grounded = proj(&heads[tipIdx*3], upAxis) < midU;
        chains.push_back(std::move(lc));
    }

    // Canonical ordering: grounded first, then by tip Z (forward) desc,
    // then by tip X asc (left-before-right under the +Y/+Z/+X convention).
    std::sort(chains.begin(), chains.end(),
        [&heads](const LegChain& a, const LegChain& b) {
            if (a.grounded != b.grounded) return a.grounded && !b.grounded;
            float za = heads[a.tipBone*3+2];
            float zb = heads[b.tipBone*3+2];
            if (std::fabs(za - zb) > 1e-4f) return za > zb;
            return heads[a.tipBone*3+0] < heads[b.tipBone*3+0];
        });

    return chains;
}

// ---- defaultGait ----------------------------------------------------------

GaitPattern defaultGait(int n) {
    GaitPattern g;
    g.dutyFactor = 0.6f;
    switch (n) {
    case 2:
        g.name = "biped_walk";
        g.phases = { 0.0f, 0.5f };
        break;
    case 4:
        // Ordered FL, FR, HL, HR. Trot: diagonal pairs FL+HR, FR+HL.
        g.name = "quadruped_trot";
        g.phases = { 0.0f, 0.5f, 0.5f, 0.0f };
        break;
    case 6:
        // Ordered FL, FR, ML, MR, RL, RR. Tripod: FL+MR+RL vs FR+ML+RR.
        g.name = "hexapod_tripod";
        g.phases = { 0.0f, 0.5f, 0.5f, 0.0f, 0.0f, 0.5f };
        break;
    case 8:
        // Ordered aL, aR, bL, bR, cL, cR, dL, dR (front to back).
        // L side: 0, 0.25, 0.5, 0.75; R side anti-phase offset by 0.5.
        g.name = "octopod_ripple";
        g.phases = { 0.00f, 0.50f, 0.25f, 0.75f, 0.50f, 0.00f, 0.75f, 0.25f };
        break;
    default:
        break; // unsupported — empty gait signals failure
    }
    return g;
}

// ---- generateLocomotionCycle ---------------------------------------------

Animation generateLocomotionCycle(const Skeleton& sk, const RigSpec& spec,
                                  const LocomotionParams& paramsIn) {
    Animation anim;
    anim.name = "locomotion";

    auto chains = identifyLegChains(sk, spec);
    std::vector<const LegChain*> grounded;
    for (const auto& c : chains) if (c.grounded) grounded.push_back(&c);
    if (grounded.empty()) return anim;

    GaitPattern gait = paramsIn.gait;
    if (gait.phases.empty()) {
        gait = defaultGait((int)grounded.size());
    }
    if (gait.phases.size() != grounded.size()) return anim;
    float duty = gait.dutyFactor;
    if (!(duty > 0.0f) || !(duty < 1.0f)) duty = 0.6f;

    float forward[3] = { paramsIn.forwardAxis[0], paramsIn.forwardAxis[1], paramsIn.forwardAxis[2] };
    float up[3]      = { paramsIn.upAxis[0],      paramsIn.upAxis[1],      paramsIn.upAxis[2] };
    vnorm(forward);
    vnorm(up);

    // Bind-pose tip world positions — the "footBase" each leg oscillates around.
    Pose bind = bindPose(sk);
    std::vector<float> bindHeads;
    allBoneHeadsWorld(sk, bind, bindHeads);
    std::vector<std::array<float, 3>> footBase(grounded.size());
    for (size_t i = 0; i < grounded.size(); ++i) {
        int tip = grounded[i]->tipBone;
        footBase[i] = { bindHeads[tip*3+0], bindHeads[tip*3+1], bindHeads[tip*3+2] };
    }

    const int nFrames = std::max(2, paramsIn.keyframesPerCycle);
    const float stride = paramsIn.strideLength;
    const float lift   = paramsIn.footLiftHeight;
    const float T      = paramsIn.cycleDuration;

    int rootBone = -1;
    for (size_t i = 0; i < sk.bones.size(); ++i) {
        if (sk.bones[i].parent < 0) { rootBone = (int)i; break; }
    }

    struct BoneKeys {
        bool tDirty = false;
        bool rDirty = false;
        std::vector<float> tVals;  // 3 per frame
        std::vector<float> rVals;  // 4 per frame
    };
    std::vector<BoneKeys> boneKeys(sk.bones.size());

    for (int f = 0; f <= nFrames; ++f) {
        const float t = (float)f / (float)nFrames; // [0, 1]

        Pose pose = bind; // fresh copy each frame

        // Subtle root-bone bob along up axis (2 bobs per cycle — one per footfall
        // in a biped walk). Simulates hip vertical oscillation.
        if (rootBone >= 0 && paramsIn.bodyBobAmplitude > 0.0f) {
            float off = paramsIn.bodyBobAmplitude * std::sin(2.0f * kPi * 2.0f * t);
            float* Tloc = &pose.data[rootBone * 10];
            Tloc[0] += up[0] * off;
            Tloc[1] += up[1] * off;
            Tloc[2] += up[2] * off;
        }

        // Per-leg foot trajectory + IK.
        for (size_t i = 0; i < grounded.size(); ++i) {
            const LegChain& lc = *grounded[i];
            float phase = t + gait.phases[i];
            phase -= std::floor(phase); // wrap to [0, 1)

            float offFwd = 0.0f;
            float offUp  = 0.0f;
            if (phase < duty) {
                // Plant: slides from +stride/2 at phase=0 to -stride/2 at phase=duty.
                float s = phase / duty;
                offFwd = stride * (0.5f - s);
            } else {
                // Swing: arcs from -stride/2 back to +stride/2 via sine lift.
                float s = (phase - duty) / (1.0f - duty);
                offFwd = stride * (-0.5f + s);
                offUp  = lift  * std::sin(kPi * s);
            }

            float target[3] = {
                footBase[i][0] + forward[0]*offFwd + up[0]*offUp,
                footBase[i][1] + forward[1]*offFwd + up[1]*offUp,
                footBase[i][2] + forward[2]*offFwd + up[2]*offUp,
            };

            if (lc.bones.size() == 3) {
                solveTwoBoneIK(sk, pose, lc.bones[0], lc.bones[1], lc.bones[2],
                               target, nullptr);
            } else if (lc.bones.size() >= 2) {
                solveFABRIK(sk, pose, lc.bones, target);
            }
        }

        // Sample all bones into per-bone buffers. Dirty-flag any bone whose
        // TRS differs from bind at this frame; only dirty bones get channels.
        for (size_t b = 0; b < sk.bones.size(); ++b) {
            const float* d = &pose.data[b * 10];
            const Bone& bone = sk.bones[b];

            auto diff3 = [](const float a[3], const float c[3]) {
                return std::fabs(a[0]-c[0]) > 1e-6f ||
                       std::fabs(a[1]-c[1]) > 1e-6f ||
                       std::fabs(a[2]-c[2]) > 1e-6f;
            };
            auto diff4 = [](const float a[4], const float c[4]) {
                return std::fabs(a[0]-c[0]) > 1e-6f ||
                       std::fabs(a[1]-c[1]) > 1e-6f ||
                       std::fabs(a[2]-c[2]) > 1e-6f ||
                       std::fabs(a[3]-c[3]) > 1e-6f;
            };
            if (diff3(&d[0], bone.localT)) boneKeys[b].tDirty = true;
            if (diff4(&d[3], bone.localR)) boneKeys[b].rDirty = true;

            boneKeys[b].tVals.push_back(d[0]);
            boneKeys[b].tVals.push_back(d[1]);
            boneKeys[b].tVals.push_back(d[2]);
            boneKeys[b].rVals.push_back(d[3]);
            boneKeys[b].rVals.push_back(d[4]);
            boneKeys[b].rVals.push_back(d[5]);
            boneKeys[b].rVals.push_back(d[6]);
        }
    }

    // Shared time track.
    std::vector<float> times(nFrames + 1);
    for (int f = 0; f <= nFrames; ++f) {
        times[f] = (float)f / (float)nFrames * T;
    }

    for (size_t b = 0; b < sk.bones.size(); ++b) {
        if (boneKeys[b].tDirty) {
            AnimChannel ch;
            ch.boneIndex = (int)b;
            ch.path = AnimChannel::Path::Translation;
            ch.interp = AnimChannel::Interp::Linear;
            ch.times = times;
            ch.values = std::move(boneKeys[b].tVals);
            anim.channels.push_back(std::move(ch));
        }
        if (boneKeys[b].rDirty) {
            AnimChannel ch;
            ch.boneIndex = (int)b;
            ch.path = AnimChannel::Path::Rotation;
            ch.interp = AnimChannel::Interp::Linear;
            ch.times = times;
            ch.values = std::move(boneKeys[b].rVals);
            anim.channels.push_back(std::move(ch));
        }
    }

    anim.duration = T;
    return anim;
}

} // namespace bromesh
