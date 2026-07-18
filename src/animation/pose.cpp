#include "bromesh/animation/pose.h"

#include <bromath/bromath.h>

#include <algorithm>
#include <cmath>
#include <cstring>

namespace bromesh {

using namespace bromath;

// ---- quaternion + matrix helpers (thin wrappers over bromath) -------------

static void quatNormalize(float* q) {
    Quat r = qnorm(Quat{q[0], q[1], q[2], q[3]});
    q[0] = r.x; q[1] = r.y; q[2] = r.z; q[3] = r.w;
}

static void quatSlerp(const float* a, const float* b, float t, float* out) {
    Quat r = qslerp(Quat{a[0], a[1], a[2], a[3]},
                    Quat{b[0], b[1], b[2], b[3]}, t);
    out[0] = r.x; out[1] = r.y; out[2] = r.z; out[3] = r.w;
}

// out = a * b (column-major 4x4)
static void matMul(const float* a, const float* b, float* out) {
    Mat4 ma, mb;
    std::memcpy(ma.data, a, 16 * sizeof(float));
    std::memcpy(mb.data, b, 16 * sizeof(float));
    Mat4 r = mmul(ma, mb);
    std::memcpy(out, r.data, 16 * sizeof(float));
}

static void composeTRS(const float* t, const float* r, const float* s, float* m) {
    Mat4 r4 = mfromTRS(Vec3{t[0], t[1], t[2]},
                       Quat{r[0], r[1], r[2], r[3]},
                       Vec3{s[0], s[1], s[2]});
    std::memcpy(m, r4.data, 16 * sizeof(float));
}

// ---- Pose construction ----------------------------------------------------

Pose bindPose(const Skeleton& skeleton) {
    Pose p;
    p.data.resize(skeleton.bones.size() * 10);
    for (size_t i = 0; i < skeleton.bones.size(); ++i) {
        float* d = &p.data[i * 10];
        const auto& b = skeleton.bones[i];
        d[0] = b.localT[0]; d[1] = b.localT[1]; d[2] = b.localT[2];
        d[3] = b.localR[0]; d[4] = b.localR[1]; d[5] = b.localR[2]; d[6] = b.localR[3];
        d[7] = b.localS[0]; d[8] = b.localS[1]; d[9] = b.localS[2];
    }
    return p;
}

// ---- Animation evaluation -------------------------------------------------

static void sampleChannel(const AnimChannel& ch, float t, float* out) {
    if (ch.times.empty()) return;
    int stride = (ch.path == AnimChannel::Path::Rotation) ? 4 : 3;

    if (t <= ch.times.front()) {
        const float* src = &ch.values[0];
        // For CUBICSPLINE the first triple is (inTangent, value, outTangent); use value.
        if (ch.interp == AnimChannel::Interp::CubicSpline) src = &ch.values[stride];
        std::memcpy(out, src, stride * sizeof(float));
        return;
    }
    if (t >= ch.times.back()) {
        size_t last = ch.times.size() - 1;
        size_t base = last * stride;
        if (ch.interp == AnimChannel::Interp::CubicSpline) base = last * stride * 3 + stride;
        std::memcpy(out, &ch.values[base], stride * sizeof(float));
        return;
    }

    // Find interval
    size_t k = 0;
    for (size_t i = 1; i < ch.times.size(); ++i) {
        if (t < ch.times[i]) { k = i - 1; break; }
    }
    float t0 = ch.times[k], t1 = ch.times[k + 1];
    float u = (t1 > t0) ? (t - t0) / (t1 - t0) : 0.0f;

    if (ch.interp == AnimChannel::Interp::Step) {
        size_t base = k * stride;
        std::memcpy(out, &ch.values[base], stride * sizeof(float));
        return;
    }

    if (ch.interp == AnimChannel::Interp::CubicSpline) {
        // Packing: (inTangent0, value0, outTangent0, inTangent1, value1, outTangent1, ...)
        size_t base0 = k * stride * 3;
        size_t base1 = (k + 1) * stride * 3;
        const float* v0 = &ch.values[base0 + stride];
        const float* m0 = &ch.values[base0 + 2 * stride];   // outTangent0
        const float* m1 = &ch.values[base1];                 // inTangent1
        const float* v1 = &ch.values[base1 + stride];
        float dt = t1 - t0;
        // bromath::chermite operates on Vec3; for stride==3 we use it directly.
        // For stride==4 (quaternion) the same Hermite polynomial applies
        // component-wise — replicate via two chermite calls on xyz then handle w.
        // Keep the expanded form so the four-component case is straightforward.
        float u2 = u*u, u3 = u2*u;
        float h00 = 2*u3 - 3*u2 + 1;
        float h10 = u3 - 2*u2 + u;
        float h01 = -2*u3 + 3*u2;
        float h11 = u3 - u2;
        for (int i = 0; i < stride; ++i)
            out[i] = h00*v0[i] + h10*dt*m0[i] + h01*v1[i] + h11*dt*m1[i];
        if (stride == 4) quatNormalize(out);
        return;
    }

    // Linear
    size_t base0 = k * stride;
    size_t base1 = (k + 1) * stride;
    if (stride == 4) {
        quatSlerp(&ch.values[base0], &ch.values[base1], u, out);
    } else {
        for (int i = 0; i < stride; ++i)
            out[i] = ch.values[base0 + i] * (1 - u) + ch.values[base1 + i] * u;
    }
}

void evaluateAnimationInto(const Skeleton& skeleton,
                           const Animation& anim,
                           float t,
                           bool loop,
                           Pose& pose) {
    if (pose.data.size() != skeleton.bones.size() * 10) {
        pose = bindPose(skeleton);
    } else {
        // Reset to bind
        for (size_t i = 0; i < skeleton.bones.size(); ++i) {
            float* d = &pose.data[i * 10];
            const auto& b = skeleton.bones[i];
            d[0]=b.localT[0]; d[1]=b.localT[1]; d[2]=b.localT[2];
            d[3]=b.localR[0]; d[4]=b.localR[1]; d[5]=b.localR[2]; d[6]=b.localR[3];
            d[7]=b.localS[0]; d[8]=b.localS[1]; d[9]=b.localS[2];
        }
    }

    float ct = t;
    if (anim.duration > 0) {
        if (loop) {
            ct = std::fmod(ct, anim.duration);
            if (ct < 0) ct += anim.duration;
        } else {
            ct = std::max(0.0f, std::min(ct, anim.duration));
        }
    }

    for (const auto& ch : anim.channels) {
        if (ch.boneIndex < 0 || ch.boneIndex >= (int)skeleton.bones.size()) continue;
        float* d = &pose.data[ch.boneIndex * 10];
        switch (ch.path) {
            case AnimChannel::Path::Translation: sampleChannel(ch, ct, d + 0); break;
            case AnimChannel::Path::Rotation:    sampleChannel(ch, ct, d + 3); break;
            case AnimChannel::Path::Scale:       sampleChannel(ch, ct, d + 7); break;
        }
    }
}

Pose evaluateAnimation(const Skeleton& skeleton, const Animation& anim, float t, bool loop) {
    Pose p = bindPose(skeleton);
    evaluateAnimationInto(skeleton, anim, t, loop, p);
    return p;
}

// ---- Blending -------------------------------------------------------------

void blendPoses(Pose& a, const Pose& b, float weight, const uint8_t* boneMask) {
    if (a.data.size() != b.data.size()) return;
    size_t bones = a.boneCount();
    for (size_t i = 0; i < bones; ++i) {
        if (boneMask && !boneMask[i]) continue;
        float* da = &a.data[i * 10];
        const float* db = &b.data[i * 10];
        // Translation
        da[0] = da[0] * (1 - weight) + db[0] * weight;
        da[1] = da[1] * (1 - weight) + db[1] * weight;
        da[2] = da[2] * (1 - weight) + db[2] * weight;
        // Rotation (slerp)
        float r[4];
        quatSlerp(&da[3], &db[3], weight, r);
        da[3] = r[0]; da[4] = r[1]; da[5] = r[2]; da[6] = r[3];
        // Scale
        da[7] = da[7] * (1 - weight) + db[7] * weight;
        da[8] = da[8] * (1 - weight) + db[8] * weight;
        da[9] = da[9] * (1 - weight) + db[9] * weight;
    }
}

void blendPosesN(const Pose* const* poses, const float* weights, size_t count,
                 Pose& out, const uint8_t* boneMask) {
    if (count == 0) return;

    // Validate shared bone count.
    size_t stride = poses[0]->data.size();
    for (size_t i = 1; i < count; ++i)
        if (poses[i]->data.size() != stride) return;

    // Normalize weights (clamped at 0); ~zero sum → first pose wins.
    float sum = 0.0f;
    for (size_t i = 0; i < count; ++i) sum += std::max(weights[i], 0.0f);

    // Reference = highest-weight source (hemisphere anchor + fill value).
    size_t ref = 0;
    for (size_t i = 1; i < count; ++i)
        if (weights[i] > weights[ref]) ref = i;

    if (out.data.size() != stride) out.data = poses[ref]->data;
    if (sum <= 1e-8f) {
        if (!boneMask) out.data = poses[0]->data;
        else {
            size_t bones = stride / 10;
            for (size_t b = 0; b < bones; ++b)
                if (boneMask[b])
                    std::copy_n(&poses[0]->data[b * 10], 10, &out.data[b * 10]);
        }
        return;
    }

    if (count == 2) {
        // Exact two-way path: identical math to blendPoses (slerp rotations)
        // so N=2 blend spaces and the two-pose crossfade agree bit-for-bit.
        float w = std::max(weights[1], 0.0f) / sum;
        size_t bones = stride / 10;
        for (size_t b = 0; b < bones; ++b) {
            if (boneMask && !boneMask[b]) continue;
            const float* d0 = &poses[0]->data[b * 10];
            const float* d1 = &poses[1]->data[b * 10];
            float* o = &out.data[b * 10];
            for (int k = 0; k < 3; ++k) o[k] = d0[k] * (1 - w) + d1[k] * w;
            quatSlerp(&d0[3], &d1[3], w, &o[3]);
            for (int k = 7; k < 10; ++k) o[k] = d0[k] * (1 - w) + d1[k] * w;
        }
        return;
    }

    size_t bones = stride / 10;
    for (size_t b = 0; b < bones; ++b) {
        if (boneMask && !boneMask[b]) continue;
        float* o = &out.data[b * 10];
        float t[3] = {0, 0, 0}, q[4] = {0, 0, 0, 0}, s[3] = {0, 0, 0};
        const float* qr = &poses[ref]->data[b * 10] + 3;
        for (size_t i = 0; i < count; ++i) {
            float w = std::max(weights[i], 0.0f) / sum;
            if (w <= 0.0f) continue;
            const float* d = &poses[i]->data[b * 10];
            t[0] += w * d[0]; t[1] += w * d[1]; t[2] += w * d[2];
            // Hemisphere-align against the reference, then weighted nlerp.
            float dot = d[3]*qr[0] + d[4]*qr[1] + d[5]*qr[2] + d[6]*qr[3];
            float sign = (dot < 0.0f) ? -w : w;
            q[0] += sign * d[3]; q[1] += sign * d[4];
            q[2] += sign * d[5]; q[3] += sign * d[6];
            s[0] += w * d[7]; s[1] += w * d[8]; s[2] += w * d[9];
        }
        o[0] = t[0]; o[1] = t[1]; o[2] = t[2];
        float len2 = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
        if (len2 > 1e-12f) {
            o[3] = q[0]; o[4] = q[1]; o[5] = q[2]; o[6] = q[3];
            quatNormalize(&o[3]);
        } else {
            // Degenerate accumulation (antipodal cancellation) — keep the
            // reference rotation rather than emitting NaN.
            o[3] = qr[0]; o[4] = qr[1]; o[5] = qr[2]; o[6] = qr[3];
        }
        o[7] = s[0]; o[8] = s[1]; o[9] = s[2];
    }
}

// ---- World + skinning matrices -------------------------------------------

void computeWorldMatrices(const Skeleton& skeleton, const Pose& pose,
                          std::vector<float>& outWorld) {
    size_t bones = skeleton.bones.size();
    outWorld.assign(bones * 16, 0.0f);

    // Identity check on rootTransform — skip the extra matMul in the common
    // case (procedurally-built skeletons, saved-out rigs with no ancestry).
    bool rootIsIdentity = true;
    for (int k = 0; k < 16 && rootIsIdentity; ++k) {
        float want = (k == 0 || k == 5 || k == 10 || k == 15) ? 1.0f : 0.0f;
        if (skeleton.rootTransform[k] != want) rootIsIdentity = false;
    }

    for (size_t i = 0; i < bones; ++i) {
        const float* d = &pose.data[i * 10];
        float local[16];
        composeTRS(&d[0], &d[3], &d[7], local);
        int parent = skeleton.bones[i].parent;
        if (parent < 0) {
            if (rootIsIdentity) {
                std::memcpy(&outWorld[i * 16], local, 16 * sizeof(float));
            } else {
                matMul(skeleton.rootTransform, local, &outWorld[i * 16]);
            }
        } else {
            matMul(&outWorld[parent * 16], local, &outWorld[i * 16]);
        }
    }
}

void computeSkinningMatrices(const Skeleton& skeleton, const Pose& pose,
                             std::vector<float>& outSkinning) {
    std::vector<float> world;
    computeWorldMatrices(skeleton, pose, world);
    size_t bones = skeleton.bones.size();
    outSkinning.assign(bones * 16, 0.0f);
    for (size_t i = 0; i < bones; ++i) {
        matMul(&world[i * 16], skeleton.bones[i].inverseBind, &outSkinning[i * 16]);
    }
}

bool socketWorldMatrix(const Skeleton& skeleton, const Pose& pose,
                       const std::string& socketName, float* outMatrix) {
    int si = skeleton.findSocket(socketName);
    if (si < 0) return false;
    const Socket& s = skeleton.sockets[si];
    if (s.bone < 0 || s.bone >= (int)skeleton.bones.size()) return false;

    std::vector<float> world;
    computeWorldMatrices(skeleton, pose, world);
    matMul(&world[s.bone * 16], s.offset, outMatrix);
    return true;
}

} // namespace bromesh
