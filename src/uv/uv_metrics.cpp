#include "bromesh/uv/uv_metrics.h"

#include <algorithm>
#include <cmath>

namespace bromesh {

// Compute area of a 3D triangle
static float triangle3DArea(const float* p0, const float* p1, const float* p2) {
    float e1[3] = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
    float e2[3] = {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};
    float cx = e1[1]*e2[2] - e1[2]*e2[1];
    float cy = e1[2]*e2[0] - e1[0]*e2[2];
    float cz = e1[0]*e2[1] - e1[1]*e2[0];
    return 0.5f * std::sqrt(cx*cx + cy*cy + cz*cz);
}

// Compute area of a 2D triangle (UV space)
static float triangle2DArea(const float* uv0, const float* uv1, const float* uv2) {
    return 0.5f * std::fabs(
        (uv1[0]-uv0[0]) * (uv2[1]-uv0[1]) -
        (uv2[0]-uv0[0]) * (uv1[1]-uv0[1])
    );
}

std::vector<UVDistortion> computeUVDistortion(const MeshData& mesh) {
    std::vector<UVDistortion> result;
    if (mesh.empty() || !mesh.hasUVs() || mesh.indices.empty()) return result;

    const size_t triCount = mesh.triangleCount();
    result.resize(triCount);

    for (size_t t = 0; t < triCount; ++t) {
        uint32_t i0 = mesh.indices[t*3+0];
        uint32_t i1 = mesh.indices[t*3+1];
        uint32_t i2 = mesh.indices[t*3+2];

        const float* p0 = &mesh.positions[i0*3];
        const float* p1 = &mesh.positions[i1*3];
        const float* p2 = &mesh.positions[i2*3];
        const float* uv0 = &mesh.uvs[i0*2];
        const float* uv1 = &mesh.uvs[i1*2];
        const float* uv2 = &mesh.uvs[i2*2];

        float area3D = triangle3DArea(p0, p1, p2);
        float area2D = triangle2DArea(uv0, uv1, uv2);

        // Area distortion: ratio of UV area to 3D area
        if (area3D > 1e-10f) {
            result[t].areaDistortion = area2D / area3D;
        }

        // L2 stretch metric (Sander et al. 2001)
        // Compute the Jacobian of the 3D->2D mapping
        // Using the singular values of the Jacobian
        if (area2D < 1e-10f || area3D < 1e-10f) {
            result[t].stretch = 1e6f; // degenerate
            result[t].angleDistortion = 1e6f;
            continue;
        }

        // 3D edges
        float e1_3d[3] = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
        float e2_3d[3] = {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};

        // UV edges
        float e1_uv[2] = {uv1[0]-uv0[0], uv1[1]-uv0[1]};
        float e2_uv[2] = {uv2[0]-uv0[0], uv2[1]-uv0[1]};

        // Inverse of the UV triangle matrix
        float det = e1_uv[0]*e2_uv[1] - e1_uv[1]*e2_uv[0];
        if (std::fabs(det) < 1e-10f) {
            result[t].stretch = 1e6f;
            result[t].angleDistortion = 1e6f;
            continue;
        }
        float invDet = 1.0f / det;

        // Jacobian columns: Ss = dp/ds, St = dp/dt
        // where s,t are UV coordinates
        float Ss[3], St[3];
        for (int c = 0; c < 3; ++c) {
            Ss[c] = invDet * ( e2_uv[1]*e1_3d[c] - e1_uv[1]*e2_3d[c]);
            St[c] = invDet * (-e2_uv[0]*e1_3d[c] + e1_uv[0]*e2_3d[c]);
        }

        // Compute coefficients of the first fundamental form
        float a = Ss[0]*Ss[0] + Ss[1]*Ss[1] + Ss[2]*Ss[2]; // dot(Ss,Ss)
        float b = Ss[0]*St[0] + Ss[1]*St[1] + Ss[2]*St[2]; // dot(Ss,St)
        float c_val = St[0]*St[0] + St[1]*St[1] + St[2]*St[2]; // dot(St,St)

        // Singular values: sigma = sqrt(eigenvalues of [[a,b],[b,c]])
        float sum = a + c_val;
        float diff = a - c_val;
        float disc = std::sqrt(std::max(0.0f, diff*diff + 4.0f*b*b));
        float sigma1_sq = 0.5f * (sum + disc); // larger
        float sigma2_sq = 0.5f * (sum - disc); // smaller
        sigma2_sq = std::max(0.0f, sigma2_sq);

        float sigma1 = std::sqrt(sigma1_sq);
        float sigma2 = std::sqrt(sigma2_sq);

        // L2 stretch: RMS of singular values (1.0 = isometric)
        result[t].stretch = std::sqrt(0.5f * (sigma1_sq + sigma2_sq));

        // Conformal (angle) distortion: ratio of singular values
        // Perfect conformal mapping has sigma1 == sigma2
        if (sigma2 > 1e-10f) {
            float ratio = sigma1 / sigma2;
            result[t].angleDistortion = ratio + 1.0f/ratio - 2.0f; // 0 = perfect
        } else {
            result[t].angleDistortion = 1e6f;
        }
    }

    return result;
}

UVMetrics measureUVQuality(const MeshData& mesh) {
    UVMetrics metrics;
    if (mesh.empty() || !mesh.hasUVs() || mesh.indices.empty()) return metrics;

    auto distortions = computeUVDistortion(mesh);
    metrics.triangleCount = distortions.size();

    // Compute 3D areas for weighting
    float totalArea3D = 0.0f;
    float totalUVArea = 0.0f;
    std::vector<float> areas3D(metrics.triangleCount);

    for (size_t t = 0; t < metrics.triangleCount; ++t) {
        uint32_t i0 = mesh.indices[t*3+0];
        uint32_t i1 = mesh.indices[t*3+1];
        uint32_t i2 = mesh.indices[t*3+2];

        areas3D[t] = triangle3DArea(
            &mesh.positions[i0*3], &mesh.positions[i1*3], &mesh.positions[i2*3]);
        totalArea3D += areas3D[t];

        totalUVArea += triangle2DArea(
            &mesh.uvs[i0*2], &mesh.uvs[i1*2], &mesh.uvs[i2*2]);
    }

    // Area-weighted averages and max
    for (size_t t = 0; t < metrics.triangleCount; ++t) {
        float w = (totalArea3D > 1e-10f) ? areas3D[t] / totalArea3D : 1.0f / metrics.triangleCount;

        metrics.avgStretch += distortions[t].stretch * w;
        metrics.maxStretch = std::max(metrics.maxStretch, distortions[t].stretch);

        metrics.avgAreaDistortion += distortions[t].areaDistortion * w;
        metrics.maxAreaDistortion = std::max(metrics.maxAreaDistortion, distortions[t].areaDistortion);

        metrics.avgAngleDistortion += distortions[t].angleDistortion * w;
        metrics.maxAngleDistortion = std::max(metrics.maxAngleDistortion, distortions[t].angleDistortion);
    }

    // UV space usage: total UV area / UV bounding box area
    float uvMin[2] = {1e10f, 1e10f};
    float uvMax[2] = {-1e10f, -1e10f};
    const size_t vCount = mesh.vertexCount();
    for (size_t v = 0; v < vCount; ++v) {
        uvMin[0] = std::min(uvMin[0], mesh.uvs[v*2+0]);
        uvMin[1] = std::min(uvMin[1], mesh.uvs[v*2+1]);
        uvMax[0] = std::max(uvMax[0], mesh.uvs[v*2+0]);
        uvMax[1] = std::max(uvMax[1], mesh.uvs[v*2+1]);
    }
    float bboxArea = (uvMax[0]-uvMin[0]) * (uvMax[1]-uvMin[1]);
    if (bboxArea > 1e-10f) {
        metrics.uvSpaceUsage = totalUVArea / bboxArea;
    }

    return metrics;
}

} // namespace bromesh
