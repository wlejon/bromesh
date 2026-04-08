#include "bromesh/analysis/sample.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>

namespace bromesh {

static float triangleArea(const float* p0, const float* p1, const float* p2) {
    float e1[3] = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
    float e2[3] = {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};
    float cx = e1[1]*e2[2] - e1[2]*e2[1];
    float cy = e1[2]*e2[0] - e1[0]*e2[2];
    float cz = e1[0]*e2[1] - e1[1]*e2[0];
    return 0.5f * std::sqrt(cx*cx + cy*cy + cz*cz);
}

std::vector<float> computeTriangleAreas(const MeshData& mesh) {
    std::vector<float> areas;
    if (mesh.empty() || mesh.indices.empty()) return areas;

    const size_t triCount = mesh.triangleCount();
    areas.resize(triCount);

    for (size_t t = 0; t < triCount; ++t) {
        uint32_t i0 = mesh.indices[t*3+0];
        uint32_t i1 = mesh.indices[t*3+1];
        uint32_t i2 = mesh.indices[t*3+2];
        areas[t] = triangleArea(
            &mesh.positions[i0*3], &mesh.positions[i1*3], &mesh.positions[i2*3]);
    }

    return areas;
}

float computeSurfaceArea(const MeshData& mesh) {
    auto areas = computeTriangleAreas(mesh);
    float total = 0.0f;
    for (float a : areas) total += a;
    return total;
}

MeshData sampleSurface(const MeshData& mesh, size_t numSamples, uint32_t seed) {
    MeshData result;
    if (mesh.empty() || mesh.indices.empty() || numSamples == 0) return result;

    const size_t triCount = mesh.triangleCount();
    const bool hasNormals = mesh.hasNormals();
    const bool hasUVs = mesh.hasUVs();

    // Compute triangle areas and build CDF for area-weighted sampling
    auto areas = computeTriangleAreas(mesh);
    std::vector<float> cdf(triCount);
    cdf[0] = areas[0];
    for (size_t t = 1; t < triCount; ++t) {
        cdf[t] = cdf[t-1] + areas[t];
    }
    float totalArea = cdf.back();
    if (totalArea < 1e-10f) return result;

    // Normalize CDF
    for (size_t t = 0; t < triCount; ++t) {
        cdf[t] /= totalArea;
    }

    // Set up RNG
    std::mt19937 rng(seed == 0 ? std::random_device{}() : seed);
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    result.positions.resize(numSamples * 3);
    if (hasNormals) result.normals.resize(numSamples * 3);
    if (hasUVs) result.uvs.resize(numSamples * 2);

    for (size_t s = 0; s < numSamples; ++s) {
        // Pick a triangle weighted by area
        float r = dist(rng);
        auto it = std::lower_bound(cdf.begin(), cdf.end(), r);
        size_t tri = static_cast<size_t>(it - cdf.begin());
        if (tri >= triCount) tri = triCount - 1;

        uint32_t i0 = mesh.indices[tri*3+0];
        uint32_t i1 = mesh.indices[tri*3+1];
        uint32_t i2 = mesh.indices[tri*3+2];

        // Random barycentric coordinates (uniform on triangle)
        float u = dist(rng);
        float v = dist(rng);
        if (u + v > 1.0f) {
            u = 1.0f - u;
            v = 1.0f - v;
        }
        float w = 1.0f - u - v;

        // Interpolate position
        for (int c = 0; c < 3; ++c) {
            result.positions[s*3+c] =
                w * mesh.positions[i0*3+c] +
                u * mesh.positions[i1*3+c] +
                v * mesh.positions[i2*3+c];
        }

        // Interpolate normal
        if (hasNormals) {
            for (int c = 0; c < 3; ++c) {
                result.normals[s*3+c] =
                    w * mesh.normals[i0*3+c] +
                    u * mesh.normals[i1*3+c] +
                    v * mesh.normals[i2*3+c];
            }
            // Normalize
            float* n = &result.normals[s*3];
            float len = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
            if (len > 1e-8f) { n[0]/=len; n[1]/=len; n[2]/=len; }
        }

        // Interpolate UVs
        if (hasUVs) {
            for (int c = 0; c < 2; ++c) {
                result.uvs[s*2+c] =
                    w * mesh.uvs[i0*2+c] +
                    u * mesh.uvs[i1*2+c] +
                    v * mesh.uvs[i2*2+c];
            }
        }
    }

    return result;
}

} // namespace bromesh
