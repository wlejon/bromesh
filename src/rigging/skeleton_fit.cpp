#include "bromesh/rigging/skeleton_fit.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <sstream>
#include <unordered_map>

namespace bromesh {

namespace {

using Vec3 = std::array<float, 3>;

Vec3 operator+(Vec3 a, Vec3 b) { return {a[0]+b[0], a[1]+b[1], a[2]+b[2]}; }
Vec3 operator-(Vec3 a, Vec3 b) { return {a[0]-b[0], a[1]-b[1], a[2]-b[2]}; }
Vec3 operator*(Vec3 a, float s) { return {a[0]*s, a[1]*s, a[2]*s}; }
float dot(Vec3 a, Vec3 b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
Vec3 cross(Vec3 a, Vec3 b) {
    return { a[1]*b[2] - a[2]*b[1],
             a[2]*b[0] - a[0]*b[2],
             a[0]*b[1] - a[1]*b[0] };
}
float length(Vec3 a) { return std::sqrt(dot(a, a)); }
Vec3 normalize(Vec3 a) {
    float l = length(a);
    if (l < 1e-12f) return {0, 1, 0};
    return a * (1.0f / l);
}

// ---- tiny expression evaluator --------------------------------------------
// Grammar (strict prefix, no whitespace handling beyond single spaces):
//   landmark:NAME
//   mid:A,B
//   lerp:A,B,T
//   offset:A,DX,DY,DZ

bool startsWith(const std::string& s, const char* p) {
    size_t n = std::strlen(p);
    return s.size() >= n && std::memcmp(s.data(), p, n) == 0;
}

std::vector<std::string> splitCsv(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == ',') { out.push_back(cur); cur.clear(); }
        else          { cur.push_back(c); }
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

bool lookupLandmark(const Landmarks& lm, const std::string& name, Vec3& out) {
    auto it = lm.points.find(name);
    if (it == lm.points.end()) return false;
    out = it->second;
    return true;
}

bool evalExpr(const std::string& expr, const Landmarks& lm, Vec3& out) {
    if (expr.empty()) return false;

    if (startsWith(expr, "landmark:")) {
        return lookupLandmark(lm, expr.substr(9), out);
    }
    if (startsWith(expr, "mid:")) {
        auto parts = splitCsv(expr.substr(4));
        if (parts.size() != 2) return false;
        Vec3 a, b;
        if (!lookupLandmark(lm, parts[0], a)) return false;
        if (!lookupLandmark(lm, parts[1], b)) return false;
        out = (a + b) * 0.5f;
        return true;
    }
    if (startsWith(expr, "lerp:")) {
        auto parts = splitCsv(expr.substr(5));
        if (parts.size() != 3) return false;
        Vec3 a, b;
        if (!lookupLandmark(lm, parts[0], a)) return false;
        if (!lookupLandmark(lm, parts[1], b)) return false;
        float t = std::stof(parts[2]);
        out = a + (b - a) * t;
        return true;
    }
    if (startsWith(expr, "offset:")) {
        auto parts = splitCsv(expr.substr(7));
        if (parts.size() != 4) return false;
        Vec3 a;
        if (!lookupLandmark(lm, parts[0], a)) return false;
        float dx = std::stof(parts[1]);
        float dy = std::stof(parts[2]);
        float dz = std::stof(parts[3]);
        out = a + Vec3{dx, dy, dz};
        return true;
    }
    return false;
}

// ---- 4x4 column-major matrix helpers --------------------------------------

void matIdentity(float* m) {
    for (int i = 0; i < 16; ++i) m[i] = 0.0f;
    m[0] = m[5] = m[10] = m[15] = 1.0f;
}

void matMul(const float* a, const float* b, float* out) {
    float tmp[16];
    for (int c = 0; c < 4; ++c) {
        for (int r = 0; r < 4; ++r) {
            float s = 0.0f;
            for (int k = 0; k < 4; ++k) s += a[k * 4 + r] * b[c * 4 + k];
            tmp[c * 4 + r] = s;
        }
    }
    std::memcpy(out, tmp, sizeof(tmp));
}

bool matInverse(const float* m, float* inv) {
    // General 4x4 inverse via cofactor expansion. Good enough for
    // rigid+scale bone matrices we deal with here.
    float a[16]; std::memcpy(a, m, sizeof(a));

    inv[0]  =  a[5]*a[10]*a[15] - a[5]*a[11]*a[14] - a[9]*a[6]*a[15] + a[9]*a[7]*a[14] + a[13]*a[6]*a[11] - a[13]*a[7]*a[10];
    inv[4]  = -a[4]*a[10]*a[15] + a[4]*a[11]*a[14] + a[8]*a[6]*a[15] - a[8]*a[7]*a[14] - a[12]*a[6]*a[11] + a[12]*a[7]*a[10];
    inv[8]  =  a[4]*a[9]*a[15]  - a[4]*a[11]*a[13] - a[8]*a[5]*a[15] + a[8]*a[7]*a[13] + a[12]*a[5]*a[11] - a[12]*a[7]*a[9];
    inv[12] = -a[4]*a[9]*a[14]  + a[4]*a[10]*a[13] + a[8]*a[5]*a[14] - a[8]*a[6]*a[13] - a[12]*a[5]*a[10] + a[12]*a[6]*a[9];

    inv[1]  = -a[1]*a[10]*a[15] + a[1]*a[11]*a[14] + a[9]*a[2]*a[15] - a[9]*a[3]*a[14] - a[13]*a[2]*a[11] + a[13]*a[3]*a[10];
    inv[5]  =  a[0]*a[10]*a[15] - a[0]*a[11]*a[14] - a[8]*a[2]*a[15] + a[8]*a[3]*a[14] + a[12]*a[2]*a[11] - a[12]*a[3]*a[10];
    inv[9]  = -a[0]*a[9]*a[15]  + a[0]*a[11]*a[13] + a[8]*a[1]*a[15] - a[8]*a[3]*a[13] - a[12]*a[1]*a[11] + a[12]*a[3]*a[9];
    inv[13] =  a[0]*a[9]*a[14]  - a[0]*a[10]*a[13] - a[8]*a[1]*a[14] + a[8]*a[2]*a[13] + a[12]*a[1]*a[10] - a[12]*a[2]*a[9];

    inv[2]  =  a[1]*a[6]*a[15]  - a[1]*a[7]*a[14]  - a[5]*a[2]*a[15] + a[5]*a[3]*a[14] + a[13]*a[2]*a[7]  - a[13]*a[3]*a[6];
    inv[6]  = -a[0]*a[6]*a[15]  + a[0]*a[7]*a[14]  + a[4]*a[2]*a[15] - a[4]*a[3]*a[14] - a[12]*a[2]*a[7]  + a[12]*a[3]*a[6];
    inv[10] =  a[0]*a[5]*a[15]  - a[0]*a[7]*a[13]  - a[4]*a[1]*a[15] + a[4]*a[3]*a[13] + a[12]*a[1]*a[7]  - a[12]*a[3]*a[5];
    inv[14] = -a[0]*a[5]*a[14]  + a[0]*a[6]*a[13]  + a[4]*a[1]*a[14] - a[4]*a[2]*a[13] - a[12]*a[1]*a[6]  + a[12]*a[2]*a[5];

    inv[3]  = -a[1]*a[6]*a[11]  + a[1]*a[7]*a[10]  + a[5]*a[2]*a[11] - a[5]*a[3]*a[10] - a[9]*a[2]*a[7]   + a[9]*a[3]*a[6];
    inv[7]  =  a[0]*a[6]*a[11]  - a[0]*a[7]*a[10]  - a[4]*a[2]*a[11] + a[4]*a[3]*a[10] + a[8]*a[2]*a[7]   - a[8]*a[3]*a[6];
    inv[11] = -a[0]*a[5]*a[11]  + a[0]*a[7]*a[9]   + a[4]*a[1]*a[11] - a[4]*a[3]*a[9]  - a[8]*a[1]*a[7]   + a[8]*a[3]*a[5];
    inv[15] =  a[0]*a[5]*a[10]  - a[0]*a[6]*a[9]   - a[4]*a[1]*a[10] + a[4]*a[2]*a[9]  + a[8]*a[1]*a[6]   - a[8]*a[2]*a[5];

    float det = a[0]*inv[0] + a[1]*inv[4] + a[2]*inv[8] + a[3]*inv[12];
    if (std::fabs(det) < 1e-20f) { matIdentity(inv); return false; }
    float invDet = 1.0f / det;
    for (int i = 0; i < 16; ++i) inv[i] *= invDet;
    return true;
}

// Build a world matrix with origin at `head`, +Y axis along `tail-head`.
// Pick +X and +Z via a stable reference vector (world-Z, falling back to
// world-X if the bone is Z-aligned).
void buildBoneWorld(Vec3 head, Vec3 tail, float* m) {
    matIdentity(m);
    Vec3 y = normalize(tail - head);

    Vec3 ref = { 0, 0, 1 };
    if (std::fabs(dot(y, ref)) > 0.95f) ref = { 1, 0, 0 };
    Vec3 x = normalize(cross(y, ref));
    Vec3 z = cross(x, y);

    m[0]  = x[0]; m[1]  = x[1]; m[2]  = x[2];
    m[4]  = y[0]; m[5]  = y[1]; m[6]  = y[2];
    m[8]  = z[0]; m[9]  = z[1]; m[10] = z[2];
    m[12] = head[0]; m[13] = head[1]; m[14] = head[2];
}

// Extract quaternion (xyzw) from the 3x3 rotation portion of m (assumed
// orthonormal). Classical Shepperd's method.
void matToQuat(const float* m, float* q) {
    float r00 = m[0], r01 = m[4], r02 = m[8];
    float r10 = m[1], r11 = m[5], r12 = m[9];
    float r20 = m[2], r21 = m[6], r22 = m[10];
    float tr = r00 + r11 + r22;
    if (tr > 0.0f) {
        float s = std::sqrt(tr + 1.0f) * 2.0f;
        q[3] = 0.25f * s;
        q[0] = (r21 - r12) / s;
        q[1] = (r02 - r20) / s;
        q[2] = (r10 - r01) / s;
    } else if (r00 > r11 && r00 > r22) {
        float s = std::sqrt(1.0f + r00 - r11 - r22) * 2.0f;
        q[3] = (r21 - r12) / s;
        q[0] = 0.25f * s;
        q[1] = (r01 + r10) / s;
        q[2] = (r02 + r20) / s;
    } else if (r11 > r22) {
        float s = std::sqrt(1.0f + r11 - r00 - r22) * 2.0f;
        q[3] = (r02 - r20) / s;
        q[0] = (r01 + r10) / s;
        q[1] = 0.25f * s;
        q[2] = (r12 + r21) / s;
    } else {
        float s = std::sqrt(1.0f + r22 - r00 - r11) * 2.0f;
        q[3] = (r10 - r01) / s;
        q[0] = (r02 + r20) / s;
        q[1] = (r12 + r21) / s;
        q[2] = 0.25f * s;
    }
}

} // namespace

// ---- main entry ------------------------------------------------------------

Skeleton fitSkeleton(const RigSpec& spec,
                     const Landmarks& landmarks,
                     const MeshData& /*mesh*/) {
    Skeleton out;

    // 1. Topologically sort bones by parent reference (spec order may be
    //    arbitrary, but parents must precede children in the output).
    std::unordered_map<std::string, int> specIndex;
    specIndex.reserve(spec.bones.size());
    for (size_t i = 0; i < spec.bones.size(); ++i) specIndex[spec.bones[i].name] = (int)i;

    std::vector<int> order;
    order.reserve(spec.bones.size());
    std::vector<int> state(spec.bones.size(), 0); // 0=unseen, 1=in-progress, 2=done
    std::vector<int> stack;

    std::function<bool(int)> visit = [&](int i) -> bool {
        if (state[i] == 2) return true;
        if (state[i] == 1) return false; // cycle
        state[i] = 1;
        const auto& p = spec.bones[i].parent;
        if (!p.empty()) {
            auto it = specIndex.find(p);
            if (it != specIndex.end() && !visit(it->second)) return false;
        }
        state[i] = 2;
        order.push_back(i);
        return true;
    };
    for (size_t i = 0; i < spec.bones.size(); ++i) visit((int)i);

    // 2. Resolve head / tail expressions against landmarks.
    std::vector<Vec3> headW(spec.bones.size(), Vec3{0,0,0});
    std::vector<Vec3> tailW(spec.bones.size(), Vec3{0,0,0});
    std::vector<bool> resolved(spec.bones.size(), false);

    // First pass: resolve what we can directly from landmarks.
    for (size_t i = 0; i < spec.bones.size(); ++i) {
        const auto& b = spec.bones[i];
        Vec3 h;
        if (evalExpr(b.head, landmarks, h)) {
            headW[i] = h;
            resolved[i] = true;
        }
    }
    // Second pass: tails. Default tail = child's head if any, else head+length*Y.
    // Find first child in spec order.
    std::vector<int> firstChild(spec.bones.size(), -1);
    for (size_t i = 0; i < spec.bones.size(); ++i) {
        const auto& p = spec.bones[i].parent;
        if (p.empty()) continue;
        auto it = specIndex.find(p);
        if (it == specIndex.end()) continue;
        int parentIdx = it->second;
        if (firstChild[parentIdx] == -1) firstChild[parentIdx] = (int)i;
    }
    for (size_t i = 0; i < spec.bones.size(); ++i) {
        const auto& b = spec.bones[i];
        Vec3 t;
        if (!b.tail.empty() && evalExpr(b.tail, landmarks, t)) {
            tailW[i] = t;
        } else if (firstChild[i] >= 0 && resolved[firstChild[i]]) {
            tailW[i] = headW[firstChild[i]];
        } else {
            float L = b.length > 0 ? b.length : 0.05f;
            tailW[i] = headW[i] + Vec3{0, L, 0};
        }
    }

    // 3. Symmetry enforcement across X=0 plane, if declared. Mean-X snap
    //    for mirrored landmark pairs — applied at bone level via paired
    //    bone names (name ends in _L/_R).
    if (spec.symmetric) {
        auto snapPair = [&](Vec3& a, Vec3& b) {
            float meanAbsX = 0.5f * (std::fabs(a[0]) + std::fabs(b[0]));
            a[0] = (a[0] >= 0 ? 1.0f : -1.0f) * meanAbsX;
            b[0] = (b[0] >= 0 ? 1.0f : -1.0f) * meanAbsX;
            float meanY = 0.5f * (a[1] + b[1]);
            float meanZ = 0.5f * (a[2] + b[2]);
            a[1] = b[1] = meanY;
            a[2] = b[2] = meanZ;
        };
        for (size_t i = 0; i < spec.bones.size(); ++i) {
            const auto& nm = spec.bones[i].name;
            if (nm.size() < 2 || nm[nm.size()-2] != '_' || nm.back() != 'L') continue;
            std::string mirrorName = nm.substr(0, nm.size()-1) + "R";
            auto it = specIndex.find(mirrorName);
            if (it == specIndex.end()) continue;
            int j = it->second;
            snapPair(headW[i], headW[j]);
            snapPair(tailW[i], tailW[j]);
        }
    }

    // 4. Build global bone world matrices + local TRS relative to parent.
    out.bones.resize(spec.bones.size());
    std::vector<float> globalM(spec.bones.size() * 16);
    std::unordered_map<std::string, int> outIndex; // name -> index into out.bones

    // Process in topo order so parent global is ready when child is visited.
    // But out.bones must also be in topo order for the runtime. Use a two-pass:
    // compute into local arrays indexed by specIndex, then copy out in order.
    std::vector<Bone> built(spec.bones.size());

    for (int si : order) {
        const auto& bspec = spec.bones[si];
        float gM[16];
        buildBoneWorld(headW[si], tailW[si], gM);
        std::memcpy(&globalM[si * 16], gM, sizeof(gM));

        // Parent global
        float parentInv[16];
        matIdentity(parentInv);
        if (!bspec.parent.empty()) {
            auto it = specIndex.find(bspec.parent);
            if (it != specIndex.end()) {
                matInverse(&globalM[it->second * 16], parentInv);
            }
        }
        float localM[16];
        matMul(parentInv, gM, localM);

        Bone b;
        b.name = bspec.name;
        b.parent = -1; // filled in after we know out-indices
        b.localT[0] = localM[12];
        b.localT[1] = localM[13];
        b.localT[2] = localM[14];
        // Extract rotation from localM (upper 3x3; scale assumed identity).
        matToQuat(localM, b.localR);
        b.localS[0] = b.localS[1] = b.localS[2] = 1.0f;

        // inverseBind = inverse(global world)
        matInverse(gM, b.inverseBind);

        built[si] = std::move(b);
    }

    // Emit in topological order.
    out.bones.clear();
    out.bones.reserve(spec.bones.size());
    for (int si : order) {
        outIndex[spec.bones[si].name] = (int)out.bones.size();
        out.bones.push_back(std::move(built[si]));
    }
    // Fix up parent indices.
    for (size_t i = 0; i < out.bones.size(); ++i) {
        // Recover parent name from spec by name lookup.
        auto sit = specIndex.find(out.bones[i].name);
        if (sit == specIndex.end()) continue;
        const auto& pname = spec.bones[sit->second].parent;
        if (pname.empty()) { out.bones[i].parent = -1; continue; }
        auto oit = outIndex.find(pname);
        out.bones[i].parent = (oit == outIndex.end()) ? -1 : oit->second;
    }

    // 5. Sockets.
    for (const auto& sd : spec.sockets) {
        auto oit = outIndex.find(sd.bone);
        if (oit == outIndex.end()) continue;
        Socket sk;
        sk.name = sd.name;
        sk.bone = oit->second;
        std::memcpy(sk.offset, sd.offset, sizeof(sk.offset));
        out.sockets.push_back(std::move(sk));
    }

    return out;
}

} // namespace bromesh
