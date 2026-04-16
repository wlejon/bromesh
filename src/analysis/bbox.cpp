#include "bromesh/analysis/bbox.h"
#include <limits>
#include <unordered_map>
#include <cmath>
#include <cstdint>
#include <tuple>
#include <utility>

namespace bromesh {

BBox computeBBox(const MeshData& mesh) {
    if (mesh.positions.empty()) return {};

    BBox bbox;
    bbox.min[0] = bbox.min[1] = bbox.min[2] =  std::numeric_limits<float>::max();
    bbox.max[0] = bbox.max[1] = bbox.max[2] = -std::numeric_limits<float>::max();

    size_t vertCount = mesh.vertexCount();
    for (size_t v = 0; v < vertCount; ++v) {
        for (int k = 0; k < 3; ++k) {
            float val = mesh.positions[v * 3 + k];
            if (val < bbox.min[k]) bbox.min[k] = val;
            if (val > bbox.max[k]) bbox.max[k] = val;
        }
    }
    return bbox;
}

bool isManifold(const MeshData& mesh) {
    if (mesh.indices.empty()) return false;

    // Weld vertices by position to get canonical vertex indices. Handles
    // meshes with split vertices (box with per-face normals, sphere seam
    // duplicates, etc.). The key is a quantized int triple so that values
    // differing only in ULP — e.g. a seam vertex at theta=2pi whose
    // computed (cos(2pi), sin(2pi)) is (0.99999994, -8.7e-8) rather than
    // (1, 0) — canonicalize with their theta=0 counterpart. Using raw
    // floats + quantized hash (an earlier attempt) created hash-only
    // collisions but not actual equality, so seam duplicates never merged.
    auto qi = [](float f) -> int32_t {
        return static_cast<int32_t>(std::round(f * 10000.0f));
    };
    struct Vec3iHash {
        size_t operator()(const std::tuple<int32_t,int32_t,int32_t>& v) const {
            size_t h = std::hash<int32_t>{}(std::get<0>(v));
            h ^= std::hash<int32_t>{}(std::get<1>(v)) + 0x9e3779b9 + (h << 6) + (h >> 2);
            h ^= std::hash<int32_t>{}(std::get<2>(v)) + 0x9e3779b9 + (h << 6) + (h >> 2);
            return h;
        }
    };

    std::unordered_map<std::tuple<int32_t,int32_t,int32_t>, uint32_t, Vec3iHash> posMap;
    std::vector<uint32_t> canonicalId(mesh.vertexCount());
    uint32_t nextId = 0;

    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        auto key = std::make_tuple(
            qi(mesh.positions[v*3]),
            qi(mesh.positions[v*3+1]),
            qi(mesh.positions[v*3+2]));
        auto it = posMap.find(key);
        if (it == posMap.end()) {
            canonicalId[v] = nextId;
            posMap[key] = nextId++;
        } else {
            canonicalId[v] = it->second;
        }
    }

    // Build edge map with canonical IDs
    std::unordered_map<uint64_t, int> edgeCount;
    size_t triCount = mesh.triangleCount();

    for (size_t t = 0; t < triCount; ++t) {
        uint32_t v[3] = {
            canonicalId[mesh.indices[t * 3 + 0]],
            canonicalId[mesh.indices[t * 3 + 1]],
            canonicalId[mesh.indices[t * 3 + 2]]
        };
        for (int e = 0; e < 3; ++e) {
            uint32_t a = v[e];
            uint32_t b = v[(e + 1) % 3];
            if (a > b) std::swap(a, b);
            uint64_t key = (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
            edgeCount[key]++;
        }
    }

    for (auto& [key, count] : edgeCount) {
        if (count != 2) return false;
    }
    return true;
}

float computeVolume(const MeshData& mesh) {
    if (!isManifold(mesh)) return 0.0f;

    // Divergence theorem: sum signed tetrahedron volumes
    double sum = 0.0;
    size_t triCount = mesh.triangleCount();
    for (size_t t = 0; t < triCount; ++t) {
        uint32_t i0 = mesh.indices[t * 3 + 0];
        uint32_t i1 = mesh.indices[t * 3 + 1];
        uint32_t i2 = mesh.indices[t * 3 + 2];

        double ax = mesh.positions[i0 * 3], ay = mesh.positions[i0 * 3 + 1], az = mesh.positions[i0 * 3 + 2];
        double bx = mesh.positions[i1 * 3], by = mesh.positions[i1 * 3 + 1], bz = mesh.positions[i1 * 3 + 2];
        double cx = mesh.positions[i2 * 3], cy = mesh.positions[i2 * 3 + 1], cz = mesh.positions[i2 * 3 + 2];

        // Signed volume of tetrahedron with origin = (1/6) * dot(a, cross(b, c))
        double crossX = by * cz - bz * cy;
        double crossY = bz * cx - bx * cz;
        double crossZ = bx * cy - by * cx;
        sum += ax * crossX + ay * crossY + az * crossZ;
    }

    return static_cast<float>(std::fabs(sum) / 6.0);
}

} // namespace bromesh
