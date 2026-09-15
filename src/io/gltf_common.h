#pragma once

#include "bromesh/io/gltf.h"

#if BROMESH_HAS_GLTF
#include "tiny_gltf.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>

namespace bromesh {

inline void matIdentity(float* m) {
    for (int i = 0; i < 16; ++i) m[i] = 0.0f;
    m[0] = m[5] = m[10] = m[15] = 1.0f;
}

// Decompose a column-major 4x4 affine matrix into T (vec3), R (quat xyzw), S (vec3).
inline void decomposeTRS(const float* m, float t[3], float r[4], float s[3]) {
    t[0] = m[12]; t[1] = m[13]; t[2] = m[14];

    float cx[3] = { m[0], m[1], m[2] };
    float cy[3] = { m[4], m[5], m[6] };
    float cz[3] = { m[8], m[9], m[10] };

    s[0] = std::sqrt(cx[0]*cx[0] + cx[1]*cx[1] + cx[2]*cx[2]);
    s[1] = std::sqrt(cy[0]*cy[0] + cy[1]*cy[1] + cy[2]*cy[2]);
    s[2] = std::sqrt(cz[0]*cz[0] + cz[1]*cz[1] + cz[2]*cz[2]);

    if (s[0] == 0 || s[1] == 0 || s[2] == 0) {
        r[0] = 0; r[1] = 0; r[2] = 0; r[3] = 1;
        return;
    }

    float rot[9] = {
        cx[0]/s[0], cx[1]/s[0], cx[2]/s[0],
        cy[0]/s[1], cy[1]/s[1], cy[2]/s[1],
        cz[0]/s[2], cz[1]/s[2], cz[2]/s[2],
    };
    // rot is column-major 3x3. Convert to quaternion (xyzw).
    float trace = rot[0] + rot[4] + rot[8];
    if (trace > 0.0f) {
        float k = std::sqrt(trace + 1.0f) * 2.0f;
        r[3] = 0.25f * k;
        r[0] = (rot[5] - rot[7]) / k;
        r[1] = (rot[6] - rot[2]) / k;
        r[2] = (rot[1] - rot[3]) / k;
    } else if (rot[0] > rot[4] && rot[0] > rot[8]) {
        float k = std::sqrt(1.0f + rot[0] - rot[4] - rot[8]) * 2.0f;
        r[3] = (rot[5] - rot[7]) / k;
        r[0] = 0.25f * k;
        r[1] = (rot[3] + rot[1]) / k;
        r[2] = (rot[6] + rot[2]) / k;
    } else if (rot[4] > rot[8]) {
        float k = std::sqrt(1.0f + rot[4] - rot[0] - rot[8]) * 2.0f;
        r[3] = (rot[6] - rot[2]) / k;
        r[0] = (rot[3] + rot[1]) / k;
        r[1] = 0.25f * k;
        r[2] = (rot[7] + rot[5]) / k;
    } else {
        float k = std::sqrt(1.0f + rot[8] - rot[0] - rot[4]) * 2.0f;
        r[3] = (rot[1] - rot[3]) / k;
        r[0] = (rot[6] + rot[2]) / k;
        r[1] = (rot[7] + rot[5]) / k;
        r[2] = 0.25f * k;
    }
}

// Pointer to raw accessor data.
inline const uint8_t* accessorData(const tinygltf::Model& model, int accessorIdx) {
    if (accessorIdx < 0 || accessorIdx >= static_cast<int>(model.accessors.size())) return nullptr;
    const auto& a = model.accessors[accessorIdx];
    if (a.bufferView < 0 || a.bufferView >= static_cast<int>(model.bufferViews.size())) return nullptr;
    const auto& bv = model.bufferViews[a.bufferView];
    if (bv.buffer < 0 || bv.buffer >= static_cast<int>(model.buffers.size())) return nullptr;
    return model.buffers[bv.buffer].data.data() + bv.byteOffset + a.byteOffset;
}

// Column-major 4x4 multiply: out = a * b.
inline void matMul4(const float* a, const float* b, float* out) {
    float tmp[16];
    for (int c = 0; c < 4; ++c) {
        for (int r = 0; r < 4; ++r) {
            tmp[c*4 + r] =
                a[0*4+r]*b[c*4+0] + a[1*4+r]*b[c*4+1] +
                a[2*4+r]*b[c*4+2] + a[3*4+r]*b[c*4+3];
        }
    }
    std::memcpy(out, tmp, 16 * sizeof(float));
}

// Build a column-major 4x4 from T, R (quat xyzw), S.
inline void trsToMat(const float t[3], const float r[4], const float s[3], float* m) {
    float xx = r[0]*r[0], yy = r[1]*r[1], zz = r[2]*r[2];
    float xy = r[0]*r[1], xz = r[0]*r[2], yz = r[1]*r[2];
    float wx = r[3]*r[0], wy = r[3]*r[1], wz = r[3]*r[2];
    m[0]  = (1.0f - 2.0f*(yy + zz)) * s[0];
    m[1]  = (2.0f*(xy + wz))        * s[0];
    m[2]  = (2.0f*(xz - wy))        * s[0];
    m[3]  = 0.0f;
    m[4]  = (2.0f*(xy - wz))        * s[1];
    m[5]  = (1.0f - 2.0f*(xx + zz)) * s[1];
    m[6]  = (2.0f*(yz + wx))        * s[1];
    m[7]  = 0.0f;
    m[8]  = (2.0f*(xz + wy))        * s[2];
    m[9]  = (2.0f*(yz - wx))        * s[2];
    m[10] = (1.0f - 2.0f*(xx + yy)) * s[2];
    m[11] = 0.0f;
    m[12] = t[0]; m[13] = t[1]; m[14] = t[2]; m[15] = 1.0f;
}

// Build the local transform matrix of a glTF node (matrix, or TRS).
inline void nodeLocalMat(const tinygltf::Node& node, float* m) {
    if (!node.matrix.empty() && node.matrix.size() == 16) {
        for (int k = 0; k < 16; ++k) m[k] = static_cast<float>(node.matrix[k]);
        return;
    }
    float t[3] = {0,0,0}, r[4] = {0,0,0,1}, s[3] = {1,1,1};
    if (node.translation.size() == 3) {
        t[0] = (float)node.translation[0];
        t[1] = (float)node.translation[1];
        t[2] = (float)node.translation[2];
    }
    if (node.rotation.size() == 4) {
        r[0] = (float)node.rotation[0];
        r[1] = (float)node.rotation[1];
        r[2] = (float)node.rotation[2];
        r[3] = (float)node.rotation[3];
    }
    if (node.scale.size() == 3) {
        s[0] = (float)node.scale[0];
        s[1] = (float)node.scale[1];
        s[2] = (float)node.scale[2];
    }
    trsToMat(t, r, s, m);
}

} // namespace bromesh
#endif // BROMESH_HAS_GLTF
