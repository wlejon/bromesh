#include "bromesh/io/splat_ply.h"
#include "ply_common.h"

#include <bromath/aabb.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

namespace bromesh {

bromath::AABB3 GaussianSplatCloud::bounds() const {
    return bromath::afromPoints(positions.data(), static_cast<int>(count()), 3);
}

namespace {

using ply_detail::PlyFormat;
using ply_detail::PlyHeader;
using ply_detail::PlyProp;

// Note: 3DGS stores base color as a band-0 SH coefficient (basis constant
// C0 = 0.5*sqrt(1/pi) ≈ 0.2820948); viewers recover color as 0.5 + C0*f_dc.
// We keep the raw coefficient in `sh` and leave SH->RGB to the rasterizer, so
// no color conversion happens in this file.

inline float sigmoidf(float x) { return 1.0f / (1.0f + std::exp(-x)); }
inline float logitf(float p) {
    // Inverse of sigmoid; clamp to keep log() finite for saturated opacities.
    p = std::min(1.0f - 1e-6f, std::max(1e-6f, p));
    return std::log(p / (1.0f - p));
}

// Count contiguous f_rest_0..f_rest_{n-1} properties present in the vertex
// element. INRIA exports them flat and channel-major: 3 * ((deg+1)^2 - 1).
int countRest(const std::vector<PlyProp>& props) {
    int n = 0;
    for (;; ++n) {
        if (ply_detail::findProp(props, "f_rest_" + std::to_string(n)) < 0) break;
    }
    return n;
}

// SH degree from the f_rest count: restTotal = 3 * ((deg+1)^2 - 1).
int shDegreeFromRest(int restTotal) {
    if (restTotal <= 0) return 0;
    int coeffsPerChannel = restTotal / 3 + 1; // (deg+1)^2
    int d = static_cast<int>(std::lround(std::sqrt(static_cast<double>(coeffsPerChannel)))) - 1;
    return std::min(3, std::max(0, d));
}

} // namespace

GaussianSplatCloud loadSplatPLY(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return {};
    std::fseek(f, 0, SEEK_END);
    long fileSize = std::ftell(f);
    std::fseek(f, 0, SEEK_SET);
    std::vector<uint8_t> fileData(fileSize);
    if (fileSize <= 0 ||
        static_cast<long>(std::fread(fileData.data(), 1, fileSize, f)) != fileSize) {
        std::fclose(f);
        return {};
    }
    std::fclose(f);

    PlyHeader header;
    if (!ply_detail::parseHeader(fileData, header)) return {};
    if (header.vertexCount == 0) return {};
    if (header.format == PlyFormat::BinaryBE) return {}; // big-endian unsupported

    const auto& props = header.vertexProps;
    const int xi = ply_detail::findProp(props, "x");
    const int yi = ply_detail::findProp(props, "y");
    const int zi = ply_detail::findProp(props, "z");
    const int oi = ply_detail::findProp(props, "opacity");
    const int s0 = ply_detail::findProp(props, "scale_0");
    const int s1 = ply_detail::findProp(props, "scale_1");
    const int s2 = ply_detail::findProp(props, "scale_2");
    const int r0 = ply_detail::findProp(props, "rot_0");
    const int r1 = ply_detail::findProp(props, "rot_1");
    const int r2 = ply_detail::findProp(props, "rot_2");
    const int r3 = ply_detail::findProp(props, "rot_3");
    const int dc0 = ply_detail::findProp(props, "f_dc_0");
    const int dc1 = ply_detail::findProp(props, "f_dc_1");
    const int dc2 = ply_detail::findProp(props, "f_dc_2");

    // Require the full splat attribute set — position, scale, rotation,
    // opacity, and DC color. Without these it is an ordinary point cloud, not
    // a Gaussian splat, and the caller should use loadPLY instead.
    if (xi < 0 || yi < 0 || zi < 0 || oi < 0 || s0 < 0 || s1 < 0 || s2 < 0 ||
        r0 < 0 || r1 < 0 || r2 < 0 || r3 < 0 || dc0 < 0 || dc1 < 0 || dc2 < 0)
        return {};

    const int restTotal = countRest(props);
    const int degree = shDegreeFromRest(restTotal);
    const int coeffs = (degree + 1) * (degree + 1); // C
    const int restPerChannel = coeffs - 1;          // f_rest per channel

    // Resolve the f_rest property indices once. Layout is channel-major:
    // f_rest_{ channel*restPerChannel + (coeff-1) }.
    std::vector<int> restIdx(static_cast<size_t>(restTotal), -1);
    for (int i = 0; i < restTotal; ++i)
        restIdx[i] = ply_detail::findProp(props, "f_rest_" + std::to_string(i));

    GaussianSplatCloud cloud;
    cloud.shDegree = degree;
    const size_t n = header.vertexCount;
    cloud.positions.resize(n * 3);
    cloud.scales.resize(n * 3);
    cloud.rotations.resize(n * 4);
    cloud.opacities.resize(n);
    cloud.sh.assign(n * static_cast<size_t>(cloud.shStride()), 0.0f);
    const int stride = cloud.shStride();

    // Pull one float property for vertex v. Binary uses precomputed offsets;
    // ASCII fills a per-vertex value row (filled by the caller below).
    auto readBinary = [&](const uint8_t* vdata, int prop) -> float {
        return ply_detail::readFloatVal(props[prop].type, vdata + header.vertexOffsets[prop]);
    };

    auto store = [&](size_t v,
                     float x, float y, float z,
                     float sc0, float sc1, float sc2,
                     float q0, float q1, float q2, float q3,
                     float op,
                     const float* dc, const float* rest) {
        cloud.positions[v * 3 + 0] = x;
        cloud.positions[v * 3 + 1] = y;
        cloud.positions[v * 3 + 2] = z;
        // Activate geometry: scales are log-space, opacity is a logit.
        cloud.scales[v * 3 + 0] = std::exp(sc0);
        cloud.scales[v * 3 + 1] = std::exp(sc1);
        cloud.scales[v * 3 + 2] = std::exp(sc2);
        // Normalize the quaternion (xyzw). 3DGS stores it w-first (rot_0=w);
        // re-emit as xyzw to match the rest of the engine.
        float qw = q0, qx = q1, qy = q2, qz = q3;
        float len = std::sqrt(qx * qx + qy * qy + qz * qz + qw * qw);
        if (len < 1e-12f) { qx = qy = qz = 0.0f; qw = 1.0f; len = 1.0f; }
        cloud.rotations[v * 4 + 0] = qx / len;
        cloud.rotations[v * 4 + 1] = qy / len;
        cloud.rotations[v * 4 + 2] = qz / len;
        cloud.rotations[v * 4 + 3] = qw / len;
        cloud.opacities[v] = sigmoidf(op);
        // DC -> coefficient 0 (interleaved RGB).
        float* dst = cloud.sh.data() + v * static_cast<size_t>(stride);
        dst[0] = dc[0];
        dst[1] = dc[1];
        dst[2] = dc[2];
        // f_rest channel-major -> interleaved by coefficient.
        for (int c = 0; c < 3; ++c)
            for (int k = 0; k < restPerChannel; ++k)
                dst[(k + 1) * 3 + c] = rest[c * restPerChannel + k];
    };

    std::vector<float> dc(3), rest(static_cast<size_t>(restTotal));

    if (header.format == PlyFormat::BinaryLE) {
        if (header.vertexStride == 0) return {};
        const uint8_t* base = fileData.data() + header.dataOffset;
        const size_t need = n * static_cast<size_t>(header.vertexStride);
        if (fileData.size() - header.dataOffset < need) return {};
        for (size_t v = 0; v < n; ++v) {
            const uint8_t* vd = base + v * header.vertexStride;
            dc[0] = readBinary(vd, dc0);
            dc[1] = readBinary(vd, dc1);
            dc[2] = readBinary(vd, dc2);
            for (int i = 0; i < restTotal; ++i) rest[i] = readBinary(vd, restIdx[i]);
            store(v,
                  readBinary(vd, xi), readBinary(vd, yi), readBinary(vd, zi),
                  readBinary(vd, s0), readBinary(vd, s1), readBinary(vd, s2),
                  readBinary(vd, r0), readBinary(vd, r1), readBinary(vd, r2), readBinary(vd, r3),
                  readBinary(vd, oi), dc.data(), rest.data());
        }
    } else {
        // ASCII: one whitespace-separated value per property, row by row.
        std::string body(reinterpret_cast<const char*>(fileData.data() + header.dataOffset),
                         fileData.size() - header.dataOffset);
        std::istringstream ss(body);
        std::vector<float> row(props.size());
        for (size_t v = 0; v < n; ++v) {
            for (size_t i = 0; i < props.size(); ++i)
                if (!(ss >> row[i])) return {};
            dc[0] = row[dc0]; dc[1] = row[dc1]; dc[2] = row[dc2];
            for (int i = 0; i < restTotal; ++i) rest[i] = row[restIdx[i]];
            store(v,
                  row[xi], row[yi], row[zi],
                  row[s0], row[s1], row[s2],
                  row[r0], row[r1], row[r2], row[r3],
                  row[oi], dc.data(), rest.data());
        }
    }

    return cloud;
}

bool saveSplatPLY(const GaussianSplatCloud& cloud, const std::string& path) {
    if (cloud.empty()) return false;

    const size_t n = cloud.count();
    const int degree = std::min(3, std::max(0, cloud.shDegree));
    const int coeffs = (degree + 1) * (degree + 1);
    const int restPerChannel = coeffs - 1;
    const int restTotal = restPerChannel * 3;
    const int stride = cloud.shStride();

    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;

    std::string header = "ply\nformat binary_little_endian 1.0\n";
    header += "element vertex " + std::to_string(n) + "\n";
    header += "property float x\nproperty float y\nproperty float z\n";
    header += "property float nx\nproperty float ny\nproperty float nz\n";
    header += "property float f_dc_0\nproperty float f_dc_1\nproperty float f_dc_2\n";
    for (int i = 0; i < restTotal; ++i)
        header += "property float f_rest_" + std::to_string(i) + "\n";
    header += "property float opacity\n";
    header += "property float scale_0\nproperty float scale_1\nproperty float scale_2\n";
    header += "property float rot_0\nproperty float rot_1\nproperty float rot_2\nproperty float rot_3\n";
    header += "end_header\n";
    std::fwrite(header.c_str(), 1, header.size(), f);

    std::vector<float> rec;
    rec.reserve(static_cast<size_t>(6 + 3 + restTotal + 1 + 3 + 4));
    for (size_t v = 0; v < n; ++v) {
        rec.clear();
        rec.push_back(cloud.positions[v * 3 + 0]);
        rec.push_back(cloud.positions[v * 3 + 1]);
        rec.push_back(cloud.positions[v * 3 + 2]);
        rec.push_back(0.0f); rec.push_back(0.0f); rec.push_back(0.0f); // normals

        const float* sh = cloud.sh.data() + v * static_cast<size_t>(stride);
        rec.push_back(sh[0]); rec.push_back(sh[1]); rec.push_back(sh[2]); // f_dc
        // interleaved-by-coefficient -> channel-major f_rest
        for (int c = 0; c < 3; ++c)
            for (int k = 0; k < restPerChannel; ++k)
                rec.push_back(sh[(k + 1) * 3 + c]);

        rec.push_back(logitf(cloud.opacities[v]));
        // Linear std-dev -> log scale. Guard against non-positive scales.
        for (int a = 0; a < 3; ++a) {
            float s = cloud.scales[v * 3 + a];
            rec.push_back(std::log(std::max(1e-12f, s)));
        }
        // xyzw in memory -> w-first rot_0..3 on disk.
        rec.push_back(cloud.rotations[v * 4 + 3]);
        rec.push_back(cloud.rotations[v * 4 + 0]);
        rec.push_back(cloud.rotations[v * 4 + 1]);
        rec.push_back(cloud.rotations[v * 4 + 2]);

        std::fwrite(rec.data(), sizeof(float), rec.size(), f);
    }

    std::fclose(f);
    return true;
}

} // namespace bromesh
