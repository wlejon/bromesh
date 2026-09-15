#include "bromesh/manipulation/splat_ops.h"
#include "bromesh/analysis/sample.h"

#include <bromath/aabb.h>
#include <bromath/mat.h>
#include <bromath/quat.h>
#include <bromath/vec.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <random>
#include <vector>

namespace bromesh {

void transformSplats(GaussianSplatCloud& cloud, const float* m) {
    if (cloud.empty() || !m) return;

    // Decompose 3x3 into rotation quaternion and axis scale lengths (column norms)
    bromath::Vec3 c0{m[0], m[1], m[2]};
    bromath::Vec3 c1{m[4], m[5], m[6]};
    bromath::Vec3 c2{m[8], m[9], m[10]};

    float sx = bromath::vlen(c0);
    float sy = bromath::vlen(c1);
    float sz = bromath::vlen(c2);
    bromath::Vec3 colScale{sx, sy, sz};

    if (sx > 1e-12f) c0 = c0 / sx;
    if (sy > 1e-12f) c1 = c1 / sy;
    if (sz > 1e-12f) c2 = c2 / sz;

    // Convert orthonormal 3x3 basis to unit quaternion (Shepperd's method)
    bromath::Quat rotMatrixQuat = bromath::qidentity();
    float tr = c0.x + c1.y + c2.z;
    if (tr > 0.0f) {
        float ss = std::sqrt(tr + 1.0f) * 2.0f;
        rotMatrixQuat.w = 0.25f * ss;
        rotMatrixQuat.x = (c1.z - c2.y) / ss;
        rotMatrixQuat.y = (c2.x - c0.z) / ss;
        rotMatrixQuat.z = (c0.y - c1.x) / ss;
    } else if (c0.x > c1.y && c0.x > c2.z) {
        float ss = std::sqrt(1.0f + c0.x - c1.y - c2.z) * 2.0f;
        rotMatrixQuat.w = (c1.z - c2.y) / ss;
        rotMatrixQuat.x = 0.25f * ss;
        rotMatrixQuat.y = (c1.x + c0.y) / ss;
        rotMatrixQuat.z = (c2.x + c0.z) / ss;
    } else if (c1.y > c2.z) {
        float ss = std::sqrt(1.0f + c1.y - c0.x - c2.z) * 2.0f;
        rotMatrixQuat.w = (c2.x - c0.z) / ss;
        rotMatrixQuat.x = (c1.x + c0.y) / ss;
        rotMatrixQuat.y = 0.25f * ss;
        rotMatrixQuat.z = (c2.y + c1.z) / ss;
    } else {
        float ss = std::sqrt(1.0f + c2.z - c0.x - c1.y) * 2.0f;
        rotMatrixQuat.w = (c0.y - c1.x) / ss;
        rotMatrixQuat.x = (c2.x + c0.z) / ss;
        rotMatrixQuat.y = (c2.y + c1.z) / ss;
        rotMatrixQuat.z = 0.25f * ss;
    }
    rotMatrixQuat = bromath::qnorm(rotMatrixQuat);

    const size_t c = cloud.count();
    const bool hasRot = !cloud.rotations.empty();
    const bool hasScale = !cloud.scales.empty();

    for (size_t v = 0; v < c; ++v) {
        float x = cloud.positions[v * 3 + 0];
        float y = cloud.positions[v * 3 + 1];
        float z = cloud.positions[v * 3 + 2];

        // Column-major: pos = m * [x, y, z, 1]
        cloud.positions[v * 3 + 0] = m[0] * x + m[4] * y + m[8] * z  + m[12];
        cloud.positions[v * 3 + 1] = m[1] * x + m[5] * y + m[9] * z  + m[13];
        cloud.positions[v * 3 + 2] = m[2] * x + m[6] * y + m[10] * z + m[14];

        if (hasRot) {
            bromath::Quat r{cloud.rotations[v * 4 + 0],
                            cloud.rotations[v * 4 + 1],
                            cloud.rotations[v * 4 + 2],
                            cloud.rotations[v * 4 + 3]};
            r = bromath::qnorm(bromath::qmul(rotMatrixQuat, r));
            cloud.rotations[v * 4 + 0] = r.x;
            cloud.rotations[v * 4 + 1] = r.y;
            cloud.rotations[v * 4 + 2] = r.z;
            cloud.rotations[v * 4 + 3] = r.w;
        }

        if (hasScale) {
            cloud.scales[v * 3 + 0] *= colScale.x;
            cloud.scales[v * 3 + 1] *= colScale.y;
            cloud.scales[v * 3 + 2] *= colScale.z;
        }
    }
}

void translateSplats(GaussianSplatCloud& cloud, float dx, float dy, float dz) {
    const size_t c = cloud.count();
    for (size_t v = 0; v < c; ++v) {
        cloud.positions[v * 3 + 0] += dx;
        cloud.positions[v * 3 + 1] += dy;
        cloud.positions[v * 3 + 2] += dz;
    }
}

void scaleSplats(GaussianSplatCloud& cloud, float sx, float sy, float sz) {
    const size_t c = cloud.count();
    const bool hasScale = !cloud.scales.empty();
    const float asx = std::abs(sx);
    const float asy = std::abs(sy);
    const float asz = std::abs(sz);

    for (size_t v = 0; v < c; ++v) {
        cloud.positions[v * 3 + 0] *= sx;
        cloud.positions[v * 3 + 1] *= sy;
        cloud.positions[v * 3 + 2] *= sz;

        if (hasScale) {
            cloud.scales[v * 3 + 0] *= asx;
            cloud.scales[v * 3 + 1] *= asy;
            cloud.scales[v * 3 + 2] *= asz;
        }
    }
}

void filterSplats(GaussianSplatCloud& cloud, const SplatFilterOptions& opts) {
    if (cloud.empty()) return;

    const size_t total = cloud.count();
    const bool hasScales = !cloud.scales.empty();
    const bool hasRots = !cloud.rotations.empty();
    const bool hasOpacities = !cloud.opacities.empty();
    const bool hasSh = !cloud.sh.empty();
    const int stride = cloud.shStride();

    size_t writeIdx = 0;
    for (size_t readIdx = 0; readIdx < total; ++readIdx) {
        if (hasOpacities && cloud.opacities[readIdx] < opts.minOpacity) {
            continue;
        }
        if (opts.cropBox) {
            bromath::Vec3 p{cloud.positions[readIdx * 3 + 0],
                            cloud.positions[readIdx * 3 + 1],
                            cloud.positions[readIdx * 3 + 2]};
            if (!bromath::acontains(*opts.cropBox, p)) {
                continue;
            }
        }
        if (opts.maxScale > 0.0f && hasScales) {
            float s0 = cloud.scales[readIdx * 3 + 0];
            float s1 = cloud.scales[readIdx * 3 + 1];
            float s2 = cloud.scales[readIdx * 3 + 2];
            if (std::max({s0, s1, s2}) > opts.maxScale) {
                continue;
            }
        }

        if (writeIdx != readIdx) {
            cloud.positions[writeIdx * 3 + 0] = cloud.positions[readIdx * 3 + 0];
            cloud.positions[writeIdx * 3 + 1] = cloud.positions[readIdx * 3 + 1];
            cloud.positions[writeIdx * 3 + 2] = cloud.positions[readIdx * 3 + 2];

            if (hasScales) {
                cloud.scales[writeIdx * 3 + 0] = cloud.scales[readIdx * 3 + 0];
                cloud.scales[writeIdx * 3 + 1] = cloud.scales[readIdx * 3 + 1];
                cloud.scales[writeIdx * 3 + 2] = cloud.scales[readIdx * 3 + 2];
            }
            if (hasRots) {
                cloud.rotations[writeIdx * 4 + 0] = cloud.rotations[readIdx * 4 + 0];
                cloud.rotations[writeIdx * 4 + 1] = cloud.rotations[readIdx * 4 + 1];
                cloud.rotations[writeIdx * 4 + 2] = cloud.rotations[readIdx * 4 + 2];
                cloud.rotations[writeIdx * 4 + 3] = cloud.rotations[readIdx * 4 + 3];
            }
            if (hasOpacities) {
                cloud.opacities[writeIdx] = cloud.opacities[readIdx];
            }
            if (hasSh) {
                std::copy_n(cloud.sh.data() + readIdx * stride,
                            stride,
                            cloud.sh.data() + writeIdx * stride);
            }
        }
        ++writeIdx;
    }

    cloud.positions.resize(writeIdx * 3);
    if (hasScales) cloud.scales.resize(writeIdx * 3);
    if (hasRots) cloud.rotations.resize(writeIdx * 4);
    if (hasOpacities) cloud.opacities.resize(writeIdx);
    if (hasSh) cloud.sh.resize(writeIdx * stride);
}

GaussianSplatCloud mergeSplats(const GaussianSplatCloud* clouds, size_t count) {
    if (!clouds || count == 0) return {};

    size_t totalCount = 0;
    int maxDegree = 0;
    bool anyScales = false;
    bool anyRots = false;
    bool anyOpacities = false;
    bool anySh = false;

    for (size_t i = 0; i < count; ++i) {
        const auto& c = clouds[i];
        totalCount += c.count();
        if (c.count() > 0) {
            maxDegree = std::max(maxDegree, std::clamp(c.shDegree, 0, 3));
            if (!c.scales.empty()) anyScales = true;
            if (!c.rotations.empty()) anyRots = true;
            if (!c.opacities.empty()) anyOpacities = true;
            if (!c.sh.empty()) anySh = true;
        }
    }

    if (totalCount == 0) return {};

    GaussianSplatCloud result;
    result.shDegree = maxDegree;
    const int dstStride = result.shStride();

    result.positions.reserve(totalCount * 3);
    if (anyScales) result.scales.reserve(totalCount * 3);
    if (anyRots) result.rotations.reserve(totalCount * 4);
    if (anyOpacities) result.opacities.reserve(totalCount);
    if (anySh) result.sh.reserve(totalCount * static_cast<size_t>(dstStride));

    for (size_t i = 0; i < count; ++i) {
        const auto& c = clouds[i];
        const size_t n = c.count();
        if (n == 0) continue;

        // Positions
        result.positions.insert(result.positions.end(), c.positions.begin(), c.positions.begin() + n * 3);

        // Scales
        if (anyScales) {
            if (!c.scales.empty()) {
                result.scales.insert(result.scales.end(), c.scales.begin(), c.scales.begin() + n * 3);
            } else {
                for (size_t v = 0; v < n; ++v) {
                    result.scales.push_back(1.0f);
                    result.scales.push_back(1.0f);
                    result.scales.push_back(1.0f);
                }
            }
        }

        // Rotations
        if (anyRots) {
            if (!c.rotations.empty()) {
                result.rotations.insert(result.rotations.end(), c.rotations.begin(), c.rotations.begin() + n * 4);
            } else {
                for (size_t v = 0; v < n; ++v) {
                    result.rotations.push_back(0.0f);
                    result.rotations.push_back(0.0f);
                    result.rotations.push_back(0.0f);
                    result.rotations.push_back(1.0f);
                }
            }
        }

        // Opacities
        if (anyOpacities) {
            if (!c.opacities.empty()) {
                result.opacities.insert(result.opacities.end(), c.opacities.begin(), c.opacities.begin() + n);
            } else {
                result.opacities.insert(result.opacities.end(), n, 1.0f);
            }
        }

        // SH
        if (anySh) {
            const int srcDegree = std::clamp(c.shDegree, 0, 3);
            const int srcStride = 3 * (srcDegree + 1) * (srcDegree + 1);
            const bool hasSrcSh = !c.sh.empty();

            for (size_t v = 0; v < n; ++v) {
                if (hasSrcSh) {
                    const float* srcPtr = c.sh.data() + v * static_cast<size_t>(srcStride);
                    result.sh.insert(result.sh.end(), srcPtr, srcPtr + srcStride);
                    if (dstStride > srcStride) {
                        result.sh.insert(result.sh.end(), static_cast<size_t>(dstStride - srcStride), 0.0f);
                    }
                } else {
                    result.sh.insert(result.sh.end(), static_cast<size_t>(dstStride), 0.0f);
                }
            }
        }
    }

    return result;
}

GaussianSplatCloud mergeSplats(const std::vector<GaussianSplatCloud>& clouds) {
    return mergeSplats(clouds.data(), clouds.size());
}

GaussianSplatCloud meshToSplats(const MeshData& mesh, const MeshToSplatsOptions& opts) {
    if (mesh.empty() || mesh.indices.empty() || opts.splatCount == 0) return {};

    const size_t triCount = mesh.triangleCount();
    auto areas = computeTriangleAreas(mesh);
    if (areas.empty()) return {};

    std::vector<float> cdf(triCount);
    cdf[0] = areas[0];
    for (size_t t = 1; t < triCount; ++t) {
        cdf[t] = cdf[t - 1] + areas[t];
    }
    float totalArea = cdf.back();
    if (totalArea < 1e-10f) return {};

    for (size_t t = 0; t < triCount; ++t) {
        cdf[t] /= totalArea;
    }

    float radius = opts.splatRadius;
    if (radius <= 0.0f) {
        radius = std::sqrt(totalArea / (static_cast<float>(opts.splatCount) * 3.14159265358979323846f)) * 1.5f;
    }

    const bool hasNormals = mesh.hasNormals();
    const bool hasColors = mesh.hasColors();

    std::mt19937 rng(opts.seed == 0 ? std::random_device{}() : opts.seed);
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    GaussianSplatCloud cloud;
    cloud.shDegree = 0;
    cloud.reserve(opts.splatCount);

    static constexpr float SH_C0 = 0.28209479177387814f;

    for (size_t s = 0; s < opts.splatCount; ++s) {
        float r = dist(rng);
        auto it = std::lower_bound(cdf.begin(), cdf.end(), r);
        size_t tri = static_cast<size_t>(it - cdf.begin());
        if (tri >= triCount) tri = triCount - 1;

        uint32_t i0 = mesh.indices[tri * 3 + 0];
        uint32_t i1 = mesh.indices[tri * 3 + 1];
        uint32_t i2 = mesh.indices[tri * 3 + 2];

        float u = dist(rng);
        float v = dist(rng);
        if (u + v > 1.0f) {
            u = 1.0f - u;
            v = 1.0f - v;
        }
        float w = 1.0f - u - v;

        // Position: sampled surface point
        float px = w * mesh.positions[i0 * 3 + 0] + u * mesh.positions[i1 * 3 + 0] + v * mesh.positions[i2 * 3 + 0];
        float py = w * mesh.positions[i0 * 3 + 1] + u * mesh.positions[i1 * 3 + 1] + v * mesh.positions[i2 * 3 + 1];
        float pz = w * mesh.positions[i0 * 3 + 2] + u * mesh.positions[i1 * 3 + 2] + v * mesh.positions[i2 * 3 + 2];

        cloud.positions.push_back(px);
        cloud.positions.push_back(py);
        cloud.positions.push_back(pz);

        // Normal: interpolated vertex normal or face normal
        float nx = 0.0f, ny = 0.0f, nz = 1.0f;
        if (hasNormals) {
            nx = w * mesh.normals[i0 * 3 + 0] + u * mesh.normals[i1 * 3 + 0] + v * mesh.normals[i2 * 3 + 0];
            ny = w * mesh.normals[i0 * 3 + 1] + u * mesh.normals[i1 * 3 + 1] + v * mesh.normals[i2 * 3 + 1];
            nz = w * mesh.normals[i0 * 3 + 2] + u * mesh.normals[i1 * 3 + 2] + v * mesh.normals[i2 * 3 + 2];
            float nlen = std::sqrt(nx * nx + ny * ny + nz * nz);
            if (nlen > 1e-8f) {
                nx /= nlen; ny /= nlen; nz /= nlen;
            } else {
                float e1x = mesh.positions[i1 * 3 + 0] - mesh.positions[i0 * 3 + 0];
                float e1y = mesh.positions[i1 * 3 + 1] - mesh.positions[i0 * 3 + 1];
                float e1z = mesh.positions[i1 * 3 + 2] - mesh.positions[i0 * 3 + 2];
                float e2x = mesh.positions[i2 * 3 + 0] - mesh.positions[i0 * 3 + 0];
                float e2y = mesh.positions[i2 * 3 + 1] - mesh.positions[i0 * 3 + 1];
                float e2z = mesh.positions[i2 * 3 + 2] - mesh.positions[i0 * 3 + 2];
                nx = e1y * e2z - e1z * e2y;
                ny = e1z * e2x - e1x * e2z;
                nz = e1x * e2y - e1y * e2x;
                float fnlen = std::sqrt(nx * nx + ny * ny + nz * nz);
                if (fnlen > 1e-8f) { nx /= fnlen; ny /= fnlen; nz /= fnlen; }
                else { nx = 0.0f; ny = 0.0f; nz = 1.0f; }
            }
        } else {
            float e1x = mesh.positions[i1 * 3 + 0] - mesh.positions[i0 * 3 + 0];
            float e1y = mesh.positions[i1 * 3 + 1] - mesh.positions[i0 * 3 + 1];
            float e1z = mesh.positions[i1 * 3 + 2] - mesh.positions[i0 * 3 + 2];
            float e2x = mesh.positions[i2 * 3 + 0] - mesh.positions[i0 * 3 + 0];
            float e2y = mesh.positions[i2 * 3 + 1] - mesh.positions[i0 * 3 + 1];
            float e2z = mesh.positions[i2 * 3 + 2] - mesh.positions[i0 * 3 + 2];
            nx = e1y * e2z - e1z * e2y;
            ny = e1z * e2x - e1x * e2z;
            nz = e1x * e2y - e1y * e2x;
            float fnlen = std::sqrt(nx * nx + ny * ny + nz * nz);
            if (fnlen > 1e-8f) { nx /= fnlen; ny /= fnlen; nz /= fnlen; }
            else { nx = 0.0f; ny = 0.0f; nz = 1.0f; }
        }

        // Rotation: quaternion rotating (0, 0, 1) to normal.
        // If normal is anti-parallel (0, 0, -1), 180° rotation about X.
        bromath::Quat rot;
        float d = nz; // dot((0,0,1), (nx, ny, nz))
        if (d > 0.999999f) {
            rot = bromath::qidentity();
        } else if (d < -0.999999f) {
            rot = bromath::Quat{1.0f, 0.0f, 0.0f, 0.0f}; // 180° rotation about X
        } else {
            float cx = -ny;
            float cy = nx;
            float cz = 0.0f;
            float s = std::sqrt((1.0f + d) * 2.0f);
            float invs = 1.0f / s;
            rot = bromath::Quat{cx * invs, cy * invs, cz * invs, s * 0.5f};
            rot = bromath::qnorm(rot);
        }
        cloud.rotations.push_back(rot.x);
        cloud.rotations.push_back(rot.y);
        cloud.rotations.push_back(rot.z);
        cloud.rotations.push_back(rot.w);

        // Scale: (radius, radius, radius * 0.1f)
        cloud.scales.push_back(radius);
        cloud.scales.push_back(radius);
        cloud.scales.push_back(radius * 0.1f);

        // Opacity
        cloud.opacities.push_back(opts.opacity);

        // SH DC term: INRIA 3DGS DC formula (color - 0.5f) / 0.28209479177387814f
        float colR = 1.0f, colG = 1.0f, colB = 1.0f;
        if (hasColors) {
            colR = w * mesh.colors[i0 * 4 + 0] + u * mesh.colors[i1 * 4 + 0] + v * mesh.colors[i2 * 4 + 0];
            colG = w * mesh.colors[i0 * 4 + 1] + u * mesh.colors[i1 * 4 + 1] + v * mesh.colors[i2 * 4 + 1];
            colB = w * mesh.colors[i0 * 4 + 2] + u * mesh.colors[i1 * 4 + 2] + v * mesh.colors[i2 * 4 + 2];
        }
        cloud.sh.push_back((colR - 0.5f) / SH_C0);
        cloud.sh.push_back((colG - 0.5f) / SH_C0);
        cloud.sh.push_back((colB - 0.5f) / SH_C0);
    }

    return cloud;
}

} // namespace bromesh
