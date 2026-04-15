#include "bromesh/animation/pose.h"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace bromesh {

// ---- quaternion + matrix helpers ------------------------------------------

static void quatNormalize(float* q) {
    float len = std::sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    if (len > 1e-8f) { q[0]/=len; q[1]/=len; q[2]/=len; q[3]/=len; }
    else             { q[0]=0; q[1]=0; q[2]=0; q[3]=1; }
}

static void quatSlerp(const float* a, const float* b, float t, float* out) {
    float cosTheta = a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
    float bb[4] = { b[0], b[1], b[2], b[3] };
    if (cosTheta < 0) { bb[0]=-bb[0]; bb[1]=-bb[1]; bb[2]=-bb[2]; bb[3]=-bb[3]; cosTheta = -cosTheta; }
    if (cosTheta > 0.9995f) {
        // Linear nlerp
        for (int i = 0; i < 4; ++i) out[i] = a[i] + t * (bb[i] - a[i]);
        quatNormalize(out);
        return;
    }
    float theta = std::acos(cosTheta);
    float sinTheta = std::sin(theta);
    float w0 = std::sin((1.0f - t) * theta) / sinTheta;
    float w1 = std::sin(t * theta) / sinTheta;
    for (int i = 0; i < 4; ++i) out[i] = a[i]*w0 + bb[i]*w1;
}

static void quatToMat(const float* q, float* m) {
    float x = q[0], y = q[1], z = q[2], w = q[3];
    float xx=x*x, yy=y*y, zz=z*z;
    float xy=x*y, xz=x*z, yz=y*z;
    float wx=w*x, wy=w*y, wz=w*z;
    m[0]  = 1 - 2*(yy+zz); m[1]  = 2*(xy+wz);     m[2]  = 2*(xz-wy);     m[3]  = 0;
    m[4]  = 2*(xy-wz);     m[5]  = 1 - 2*(xx+zz); m[6]  = 2*(yz+wx);     m[7]  = 0;
    m[8]  = 2*(xz+wy);     m[9]  = 2*(yz-wx);     m[10] = 1 - 2*(xx+yy); m[11] = 0;
    m[12] = 0;             m[13] = 0;             m[14] = 0;             m[15] = 1;
}

// out = a * b (column-major 4x4)
static void matMul(const float* a, const float* b, float* out) {
    float r[16];
    for (int c = 0; c < 4; ++c) {
        for (int row = 0; row < 4; ++row) {
            float s = 0;
            for (int k = 0; k < 4; ++k)
                s += a[k * 4 + row] * b[c * 4 + k];
            r[c * 4 + row] = s;
        }
    }
    std::memcpy(out, r, 16 * sizeof(float));
}

static void composeTRS(const float* t, const float* r, const float* s, float* m) {
    float rot[16];
    quatToMat(r, rot);
    rot[0] *= s[0]; rot[1] *= s[0]; rot[2] *= s[0];
    rot[4] *= s[1]; rot[5] *= s[1]; rot[6] *= s[1];
    rot[8] *= s[2]; rot[9] *= s[2]; rot[10] *= s[2];
    rot[12] = t[0]; rot[13] = t[1]; rot[14] = t[2];
    std::memcpy(m, rot, 16 * sizeof(float));
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
        float u2 = u*u, u3 = u2*u;
        float h00 = 2*u3 - 3*u2 + 1;
        float h10 = u3 - 2*u2 + u;
        float h01 = -2*u3 + 3*u2;
        float h11 = u3 - u2;
        float dt = t1 - t0;
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

// ---- World + skinning matrices -------------------------------------------

void computeWorldMatrices(const Skeleton& skeleton, const Pose& pose,
                          std::vector<float>& outWorld) {
    size_t bones = skeleton.bones.size();
    outWorld.assign(bones * 16, 0.0f);

    for (size_t i = 0; i < bones; ++i) {
        const float* d = &pose.data[i * 10];
        float local[16];
        composeTRS(&d[0], &d[3], &d[7], local);
        int parent = skeleton.bones[i].parent;
        if (parent < 0) {
            std::memcpy(&outWorld[i * 16], local, 16 * sizeof(float));
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
