#pragma once

#include "bromesh/io/gltf.h"

#if BROMESH_HAS_GLTF
#include "tiny_gltf.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <utility>
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

// One component of a vertex attribute element as a float: honours the
// component type and `normalized` (glTF's unorm/snorm rules). `p` points at
// the component.
inline float readComponentAsFloat(const uint8_t* p, int componentType, bool normalized) {
    switch (componentType) {
        case TINYGLTF_COMPONENT_TYPE_FLOAT: {
            float f; std::memcpy(&f, p, 4); return f;
        }
        case TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE:
            return normalized ? p[0] / 255.0f : static_cast<float>(p[0]);
        case TINYGLTF_COMPONENT_TYPE_BYTE: {
            float v = static_cast<float>(static_cast<int8_t>(p[0]));
            return normalized ? std::max(v / 127.0f, -1.0f) : v;
        }
        case TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT: {
            uint16_t u; std::memcpy(&u, p, 2);
            return normalized ? u / 65535.0f : static_cast<float>(u);
        }
        case TINYGLTF_COMPONENT_TYPE_SHORT: {
            int16_t s; std::memcpy(&s, p, 2);
            return normalized ? std::max(s / 32767.0f, -1.0f) : static_cast<float>(s);
        }
        case TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT: {
            uint32_t u; std::memcpy(&u, p, 4); return static_cast<float>(u);
        }
        default: return 0.0f;
    }
}

// Byte span of one accessor element walk: base pointer, per-element stride,
// and the component size, bounds-checked against the buffer. False when the
// accessor has no data or would read past its buffer.
struct AccessorWalk {
    const uint8_t* base = nullptr;
    size_t stride = 0;
    size_t compSize = 0;
    int comps = 0;
    size_t count = 0;
};
inline bool accessorWalk(const tinygltf::Model& model, int accessorIdx, AccessorWalk& w) {
    if (accessorIdx < 0 || accessorIdx >= static_cast<int>(model.accessors.size())) return false;
    const auto& a = model.accessors[accessorIdx];
    if (a.bufferView < 0 || a.bufferView >= static_cast<int>(model.bufferViews.size())) return false;
    const auto& bv = model.bufferViews[a.bufferView];
    if (bv.buffer < 0 || bv.buffer >= static_cast<int>(model.buffers.size())) return false;
    const int compSize = tinygltf::GetComponentSizeInBytes(static_cast<uint32_t>(a.componentType));
    const int comps = tinygltf::GetNumComponentsInType(static_cast<uint32_t>(a.type));
    if (compSize <= 0 || comps <= 0) return false;
    const size_t elemSize = static_cast<size_t>(compSize) * comps;
    const size_t stride = bv.byteStride > 0 ? bv.byteStride : elemSize;
    const auto& buf = model.buffers[bv.buffer].data;
    const size_t start = bv.byteOffset + a.byteOffset;
    if (a.count > 0) {
        const size_t end = start + stride * (a.count - 1) + elemSize;
        if (end > buf.size() || end > bv.byteOffset + bv.byteLength) return false;
    }
    w.base = buf.data() + start;
    w.stride = stride;
    w.compSize = static_cast<size_t>(compSize);
    w.comps = comps;
    w.count = a.count;
    return true;
}

// A vertex attribute as `outComps` floats per element, whatever its storage:
// interleaved (byteStride) or tight, float or (normalized) integer. Missing
// source components read 0 (the 4th of a VEC3 colour is left to the caller).
inline bool readAccessorFloats(const tinygltf::Model& model, int accessorIdx,
                               int outComps, std::vector<float>& out) {
    AccessorWalk w;
    if (!accessorWalk(model, accessorIdx, w)) return false;
    const auto& a = model.accessors[accessorIdx];
    out.assign(w.count * outComps, 0.0f);
    const int n = std::min(w.comps, outComps);
    for (size_t i = 0; i < w.count; ++i) {
        const uint8_t* e = w.base + i * w.stride;
        for (int c = 0; c < n; ++c)
            out[i * outComps + c] = readComponentAsFloat(e + c * w.compSize, a.componentType, a.normalized);
    }
    return true;
}

// An integer vertex attribute (JOINTS_n) as `outComps` uint32 per element.
inline bool readAccessorUints(const tinygltf::Model& model, int accessorIdx,
                              int outComps, std::vector<uint32_t>& out) {
    AccessorWalk w;
    if (!accessorWalk(model, accessorIdx, w)) return false;
    const auto& a = model.accessors[accessorIdx];
    out.assign(w.count * outComps, 0u);
    const int n = std::min(w.comps, outComps);
    for (size_t i = 0; i < w.count; ++i) {
        const uint8_t* e = w.base + i * w.stride;
        for (int c = 0; c < n; ++c) {
            const uint8_t* p = e + c * w.compSize;
            uint32_t v = 0;
            switch (a.componentType) {
                case TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE: v = p[0]; break;
                case TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT: { uint16_t u; std::memcpy(&u, p, 2); v = u; break; }
                case TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT: std::memcpy(&v, p, 4); break;
                default: v = static_cast<uint32_t>(readComponentAsFloat(p, a.componentType, false)); break;
            }
            out[i * outComps + c] = v;
        }
    }
    return true;
}

// Fold any number of glTF influence sets (JOINTS_n / WEIGHTS_n, each 4 wide)
// into SkinData's fixed 4 influences per vertex: the 4 heaviest across all
// sets are kept and renormalized to sum to 1. A joint that appears in more
// than one set has its weights summed first. `sets` pairs a joints and a
// weights array, both 4 per vertex.
inline void mergeInfluenceSets(
        const std::vector<std::pair<std::vector<uint32_t>, std::vector<float>>>& sets,
        std::vector<uint32_t>& outJoints, std::vector<float>& outWeights) {
    outJoints.clear();
    outWeights.clear();
    if (sets.empty()) return;
    if (sets.size() == 1) {
        // Already 4 wide: keep the file's slots and values exactly.
        const size_t n = std::min(sets[0].first.size(), sets[0].second.size()) / 4 * 4;
        outJoints.assign(sets[0].first.begin(), sets[0].first.begin() + n);
        outWeights.assign(sets[0].second.begin(), sets[0].second.begin() + n);
        return;
    }
    size_t vcount = SIZE_MAX;
    for (const auto& s : sets)
        vcount = std::min(vcount, std::min(s.first.size(), s.second.size()) / 4);
    if (vcount == SIZE_MAX) return;
    outJoints.assign(vcount * 4, 0u);
    outWeights.assign(vcount * 4, 0.0f);
    std::vector<std::pair<float, uint32_t>> cand;
    cand.reserve(sets.size() * 4);
    for (size_t v = 0; v < vcount; ++v) {
        cand.clear();
        for (const auto& s : sets) {
            for (int k = 0; k < 4; ++k) {
                const float w = s.second[v * 4 + k];
                if (!(w > 0.0f)) continue;  // drops zeros, negatives and NaN
                const uint32_t j = s.first[v * 4 + k];
                bool merged = false;
                for (auto& c : cand) {
                    if (c.second == j) { c.first += w; merged = true; break; }
                }
                if (!merged) cand.emplace_back(w, j);
            }
        }
        // Heaviest first; ties keep set/slot order.
        std::stable_sort(cand.begin(), cand.end(),
                         [](const auto& a, const auto& b) { return a.first > b.first; });
        const size_t keep = std::min<size_t>(cand.size(), 4);
        float sum = 0.0f;
        for (size_t k = 0; k < keep; ++k) sum += cand[k].first;
        const float inv = sum > 0.0f ? 1.0f / sum : 0.0f;
        for (size_t k = 0; k < keep; ++k) {
            outJoints[v * 4 + k] = cand[k].second;
            outWeights[v * 4 + k] = cand[k].first * inv;
        }
    }
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
