#include "bromesh/animation/ik.h"

#include <bromath/bromath.h>

#include <algorithm>
#include <cmath>
#include <cstring>

namespace bromesh {

using namespace bromath;

// ---- float-pointer adapters over bromath ---------------------------------
// The IK routines operate on raw float buffers (consumer-friendly when
// interoperating with Pose data) — these are thin wrappers that pack into
// bromath types, call the canonical op, and write back.

static float v3len(const float* v) {
    return vlen(Vec3{v[0], v[1], v[2]});
}
static void v3sub(const float* a, const float* b, float* o) {
    o[0]=a[0]-b[0]; o[1]=a[1]-b[1]; o[2]=a[2]-b[2];
}
static void v3scale(const float* a, float s, float* o) {
    o[0]=a[0]*s; o[1]=a[1]*s; o[2]=a[2]*s;
}
static float v3dot(const float* a, const float* b) {
    return vdot(Vec3{a[0],a[1],a[2]}, Vec3{b[0],b[1],b[2]});
}
static void v3cross(const float* a, const float* b, float* o) {
    Vec3 r = vcross(Vec3{a[0],a[1],a[2]}, Vec3{b[0],b[1],b[2]});
    o[0]=r.x; o[1]=r.y; o[2]=r.z;
}
static void v3norm(float* v) {
    Vec3 r = vnorm(Vec3{v[0], v[1], v[2]});
    v[0]=r.x; v[1]=r.y; v[2]=r.z;
}

// Quaternion product p*q (xyzw) -> out
static void quatMulRaw(const float* p, const float* q, float* out) {
    Quat r = qmul(Quat{p[0],p[1],p[2],p[3]}, Quat{q[0],q[1],q[2],q[3]});
    out[0]=r.x; out[1]=r.y; out[2]=r.z; out[3]=r.w;
}

static void quatConj(const float* q, float* out) {
    Quat r = qconjugate(Quat{q[0],q[1],q[2],q[3]});
    out[0]=r.x; out[1]=r.y; out[2]=r.z; out[3]=r.w;
}

static void quatNormalize(float* q) {
    Quat r = qnorm(Quat{q[0],q[1],q[2],q[3]});
    q[0]=r.x; q[1]=r.y; q[2]=r.z; q[3]=r.w;
}

// Rotate vector v by quaternion q -> out
static void quatRotateVec(const float* q, const float* v, float* out) {
    Vec3 r = qrotate(Quat{q[0],q[1],q[2],q[3]}, Vec3{v[0],v[1],v[2]});
    out[0]=r.x; out[1]=r.y; out[2]=r.z;
}

// Quaternion that rotates `from` to `to` (both normalized vec3).
static void quatFromToRaw(const float* from, const float* to, float* out) {
    Quat r = qfromTo(Vec3{from[0],from[1],from[2]}, Vec3{to[0],to[1],to[2]});
    out[0]=r.x; out[1]=r.y; out[2]=r.z; out[3]=r.w;
}

// Gather accumulated world rotation and world position by walking bone chain
// from root, using local TRS in pose.
static void accumulateBoneWorld(const Skeleton& skeleton, const Pose& pose,
                                int bone, float outPos[3], float outRot[4]) {
    // Build list of ancestors
    int path[64];
    int depth = 0;
    int cur = bone;
    while (cur >= 0 && depth < 64) {
        path[depth++] = cur;
        cur = skeleton.bones[cur].parent;
    }
    // Walk root->bone
    outPos[0]=outPos[1]=outPos[2]=0;
    outRot[0]=0; outRot[1]=0; outRot[2]=0; outRot[3]=1;
    for (int i = depth - 1; i >= 0; --i) {
        int b = path[i];
        const float* d = &pose.data[b * 10];
        // world_pos = parent_pos + parent_rot * local_t  (ignore scale for IK sanity)
        float rotated[3];
        quatRotateVec(outRot, &d[0], rotated);
        outPos[0] += rotated[0];
        outPos[1] += rotated[1];
        outPos[2] += rotated[2];
        // world_rot = parent_rot * local_rot
        float q[4];
        quatMulRaw(outRot, &d[3], q);
        quatNormalize(q);
        outRot[0]=q[0]; outRot[1]=q[1]; outRot[2]=q[2]; outRot[3]=q[3];
    }
}

// Apply a world-space rotation `qDelta` to bone's local rotation: newLocal =
// inverse(parentWorld) * qDelta * oldWorld  (relative to parent)
// Equivalent: newLocal = inv(parentWorld) * qDelta * parentWorld * oldLocal
static void applyWorldRotationDelta(const Skeleton& skeleton, Pose& pose,
                                     int bone, const float qDelta[4]) {
    int parent = skeleton.bones[bone].parent;
    float pPos[3], pRot[4] = {0,0,0,1};
    if (parent >= 0) accumulateBoneWorld(skeleton, pose, parent, pPos, pRot);

    float pRotConj[4];
    quatConj(pRot, pRotConj);

    // tmp = pRotConj * qDelta
    float tmp[4];
    quatMulRaw(pRotConj, qDelta, tmp);
    // tmp = tmp * pRot
    quatMulRaw(tmp, pRot, tmp);
    // newLocal = tmp * oldLocal
    float* local = &pose.data[bone * 10 + 3];
    float result[4];
    quatMulRaw(tmp, local, result);
    quatNormalize(result);
    local[0]=result[0]; local[1]=result[1]; local[2]=result[2]; local[3]=result[3];
}

// ---- Two-bone IK ---------------------------------------------------------

bool solveTwoBoneIK(const Skeleton& skeleton,
                    Pose& pose,
                    int rootBone, int midBone, int endBone,
                    const float targetWorld[3],
                    const float* poleWorld) {
    if (rootBone < 0 || midBone < 0 || endBone < 0) return false;
    if ((int)skeleton.bones.size() <= std::max({rootBone, midBone, endBone})) return false;
    if (skeleton.bones[midBone].parent != rootBone) return false;
    if (skeleton.bones[endBone].parent != midBone) return false;

    // Current world positions
    float rootPos[3], rootRot[4];
    accumulateBoneWorld(skeleton, pose, rootBone, rootPos, rootRot);
    float midPos[3], midRot[4];
    accumulateBoneWorld(skeleton, pose, midBone, midPos, midRot);
    float endPos[3], endRot[4];
    accumulateBoneWorld(skeleton, pose, endBone, endPos, endRot);

    float a[3]; v3sub(midPos, rootPos, a);
    float b[3]; v3sub(endPos, midPos, b);
    float c[3]; v3sub((float*)targetWorld, rootPos, c);

    float aLen = v3len(a), bLen = v3len(b), cLen = v3len(c);
    if (aLen < 1e-6f || bLen < 1e-6f) return false;

    // Clamp target so the triangle is realizable.
    float maxReach = aLen + bLen;
    if (cLen > maxReach) { v3norm(c); v3scale(c, maxReach - 1e-4f, c); cLen = maxReach - 1e-4f; }
    float minReach = std::fabs(aLen - bLen);
    if (cLen < minReach) { v3norm(c); v3scale(c, minReach + 1e-4f, c); cLen = minReach + 1e-4f; }

    // Compute desired elbow position in the plane through root, target, and
    // either the pole (if given) or the current elbow position.
    // 1) Base vector: unit vector from root to target.
    float dHat[3] = { c[0]/cLen, c[1]/cLen, c[2]/cLen };

    // 2) Pole hint direction: from root, perpendicular to dHat, pointing toward
    //    either the pole or the current elbow.
    float hint[3];
    if (poleWorld) {
        v3sub((float*)poleWorld, rootPos, hint);
    } else {
        v3sub(midPos, rootPos, hint);
    }
    // Remove component along dHat
    float projLen = v3dot(hint, dHat);
    hint[0] -= dHat[0] * projLen;
    hint[1] -= dHat[1] * projLen;
    hint[2] -= dHat[2] * projLen;
    if (v3len(hint) < 1e-5f) {
        // Chain was fully straight; choose any perpendicular
        float up[3] = { 0, 1, 0 };
        v3cross(dHat, up, hint);
        if (v3len(hint) < 1e-4f) { up[0]=1; up[1]=0; up[2]=0; v3cross(dHat, up, hint); }
    }
    v3norm(hint);

    // 3) Distance from root along dHat to elbow (projection) via triangle math:
    //    x = (a^2 - b^2 + c^2) / (2c); y = sqrt(a^2 - x^2)
    float x = (aLen*aLen - bLen*bLen + cLen*cLen) / (2.0f * cLen);
    float y2 = aLen*aLen - x*x;
    float y = (y2 > 0) ? std::sqrt(y2) : 0.0f;

    float desiredMid[3] = {
        rootPos[0] + dHat[0]*x + hint[0]*y,
        rootPos[1] + dHat[1]*x + hint[1]*y,
        rootPos[2] + dHat[2]*x + hint[2]*y,
    };

    // --- Step A: rotate root so its child (elbow) moves to desiredMid ---
    float curRootToMid[3]; v3sub(midPos, rootPos, curRootToMid); v3norm(curRootToMid);
    float newRootToMid[3]; v3sub(desiredMid, rootPos, newRootToMid); v3norm(newRootToMid);
    float qA[4];
    quatFromToRaw(curRootToMid, newRootToMid, qA);
    applyWorldRotationDelta(skeleton, pose, rootBone, qA);

    // Refresh world positions
    accumulateBoneWorld(skeleton, pose, midBone, midPos, midRot);
    accumulateBoneWorld(skeleton, pose, endBone, endPos, endRot);

    // --- Step B: rotate mid so its child (wrist) moves to the clamped target ---
    float clampedTarget[3] = { rootPos[0] + c[0], rootPos[1] + c[1], rootPos[2] + c[2] };
    float curMidToEnd[3]; v3sub(endPos, midPos, curMidToEnd); v3norm(curMidToEnd);
    float newMidToEnd[3]; v3sub(clampedTarget, midPos, newMidToEnd); v3norm(newMidToEnd);
    float qB[4];
    quatFromToRaw(curMidToEnd, newMidToEnd, qB);
    applyWorldRotationDelta(skeleton, pose, midBone, qB);

    return true;
}

// ---- FABRIK ---------------------------------------------------------------

bool solveFABRIK(const Skeleton& skeleton,
                 Pose& pose,
                 const std::vector<int>& chain,
                 const float targetWorld[3],
                 int iterations,
                 float tolerance) {
    size_t n = chain.size();
    if (n < 2) return false;
    for (int b : chain) if (b < 0 || b >= (int)skeleton.bones.size()) return false;

    // Gather world positions and segment lengths from current pose.
    std::vector<float> pos(n * 3);
    for (size_t i = 0; i < n; ++i) {
        float p[3], r[4];
        accumulateBoneWorld(skeleton, pose, chain[i], p, r);
        pos[i*3+0] = p[0]; pos[i*3+1] = p[1]; pos[i*3+2] = p[2];
    }
    std::vector<float> segLen(n - 1);
    float totalLen = 0;
    for (size_t i = 0; i + 1 < n; ++i) {
        float d[3] = {
            pos[(i+1)*3+0] - pos[i*3+0],
            pos[(i+1)*3+1] - pos[i*3+1],
            pos[(i+1)*3+2] - pos[i*3+2],
        };
        segLen[i] = v3len(d);
        totalLen += segLen[i];
    }

    float rootPos[3] = { pos[0], pos[1], pos[2] };
    float dRoot[3]; v3sub((float*)targetWorld, rootPos, dRoot);

    // If target out of reach, stretch chain toward it
    if (v3len(dRoot) > totalLen) {
        float dir[3] = { dRoot[0], dRoot[1], dRoot[2] };
        v3norm(dir);
        for (size_t i = 1; i < n; ++i) {
            float off[3];
            v3scale(dir, segLen[i-1], off);
            pos[i*3+0] = pos[(i-1)*3+0] + off[0];
            pos[i*3+1] = pos[(i-1)*3+1] + off[1];
            pos[i*3+2] = pos[(i-1)*3+2] + off[2];
        }
    } else {
        // Iterative forward/backward
        for (int iter = 0; iter < iterations; ++iter) {
            // Backward: tip to root
            pos[(n-1)*3+0] = targetWorld[0];
            pos[(n-1)*3+1] = targetWorld[1];
            pos[(n-1)*3+2] = targetWorld[2];
            for (int i = (int)n - 2; i >= 0; --i) {
                float d[3] = {
                    pos[i*3+0] - pos[(i+1)*3+0],
                    pos[i*3+1] - pos[(i+1)*3+1],
                    pos[i*3+2] - pos[(i+1)*3+2],
                };
                v3norm(d);
                pos[i*3+0] = pos[(i+1)*3+0] + d[0] * segLen[i];
                pos[i*3+1] = pos[(i+1)*3+1] + d[1] * segLen[i];
                pos[i*3+2] = pos[(i+1)*3+2] + d[2] * segLen[i];
            }
            // Forward: root to tip
            pos[0] = rootPos[0]; pos[1] = rootPos[1]; pos[2] = rootPos[2];
            for (size_t i = 0; i + 1 < n; ++i) {
                float d[3] = {
                    pos[(i+1)*3+0] - pos[i*3+0],
                    pos[(i+1)*3+1] - pos[i*3+1],
                    pos[(i+1)*3+2] - pos[i*3+2],
                };
                v3norm(d);
                pos[(i+1)*3+0] = pos[i*3+0] + d[0] * segLen[i];
                pos[(i+1)*3+1] = pos[i*3+1] + d[1] * segLen[i];
                pos[(i+1)*3+2] = pos[i*3+2] + d[2] * segLen[i];
            }
            float tipDiff[3] = {
                pos[(n-1)*3+0] - targetWorld[0],
                pos[(n-1)*3+1] - targetWorld[1],
                pos[(n-1)*3+2] - targetWorld[2],
            };
            if (v3len(tipDiff) < tolerance) break;
        }
    }

    // Convert positions back into rotations per bone.
    // For each bone in the chain except the last: rotate so that its child's
    // original direction points toward the new child position.
    for (size_t i = 0; i + 1 < n; ++i) {
        int b = chain[i];
        // Old world direction from b to child
        float oldChild[3], oldB[3], oldChildR[4], oldBR[4];
        accumulateBoneWorld(skeleton, pose, b, oldB, oldBR);
        accumulateBoneWorld(skeleton, pose, chain[i+1], oldChild, oldChildR);
        float oldDir[3]; v3sub(oldChild, oldB, oldDir); v3norm(oldDir);

        float newDir[3] = {
            pos[(i+1)*3+0] - pos[i*3+0],
            pos[(i+1)*3+1] - pos[i*3+1],
            pos[(i+1)*3+2] - pos[i*3+2],
        };
        v3norm(newDir);

        float q[4];
        quatFromToRaw(oldDir, newDir, q);
        applyWorldRotationDelta(skeleton, pose, b, q);
    }

    return true;
}

// ---- Look-at -------------------------------------------------------------

bool solveLookAt(const Skeleton& skeleton,
                 Pose& pose,
                 int bone,
                 const float targetWorld[3],
                 const float localForward[3],
                 const float localUp[3]) {
    if (bone < 0 || bone >= (int)skeleton.bones.size()) return false;
    float fwd[3] = { 0, 0, 1 };
    float up[3]  = { 0, 1, 0 };
    (void)up;
    if (localForward) { fwd[0]=localForward[0]; fwd[1]=localForward[1]; fwd[2]=localForward[2]; v3norm(fwd); }
    if (localUp)      { up[0]=localUp[0]; up[1]=localUp[1]; up[2]=localUp[2]; v3norm(up); }

    float bonePos[3], boneRot[4];
    accumulateBoneWorld(skeleton, pose, bone, bonePos, boneRot);

    // Current world forward = boneRot * localForward
    float curFwd[3];
    quatRotateVec(boneRot, fwd, curFwd);

    float desiredFwd[3]; v3sub((float*)targetWorld, bonePos, desiredFwd);
    if (v3len(desiredFwd) < 1e-6f) return false;
    v3norm(desiredFwd);

    float q[4];
    quatFromToRaw(curFwd, desiredFwd, q);
    applyWorldRotationDelta(skeleton, pose, bone, q);
    return true;
}

} // namespace bromesh
