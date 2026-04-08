#include "bromesh/analysis/raycast.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace bromesh {

// Möller-Trumbore ray-triangle intersection returning barycentric coords
static bool rayTriangle(
    const float* orig, const float* dir,
    const float* v0, const float* v1, const float* v2,
    float& t, float& u, float& v) {

    float e1[3] = {v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]};
    float e2[3] = {v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]};

    float h[3] = {
        dir[1]*e2[2] - dir[2]*e2[1],
        dir[2]*e2[0] - dir[0]*e2[2],
        dir[0]*e2[1] - dir[1]*e2[0]
    };
    float a = e1[0]*h[0] + e1[1]*h[1] + e1[2]*h[2];
    if (std::fabs(a) < 1e-10f) return false;

    float f = 1.0f / a;
    float s[3] = {orig[0]-v0[0], orig[1]-v0[1], orig[2]-v0[2]};
    u = f * (s[0]*h[0] + s[1]*h[1] + s[2]*h[2]);
    if (u < 0.0f || u > 1.0f) return false;

    float q[3] = {
        s[1]*e1[2] - s[2]*e1[1],
        s[2]*e1[0] - s[0]*e1[2],
        s[0]*e1[1] - s[1]*e1[0]
    };
    v = f * (dir[0]*q[0] + dir[1]*q[1] + dir[2]*q[2]);
    if (v < 0.0f || u + v > 1.0f) return false;

    t = f * (e2[0]*q[0] + e2[1]*q[1] + e2[2]*q[2]);
    return t > 0.0f;
}

static void fillHit(RayHit& hit, const MeshData& mesh,
                     const float* origin, const float* dir,
                     float t, float u, float v, uint32_t triIdx) {
    hit.hit = true;
    hit.distance = t;
    hit.triangleIndex = triIdx;
    hit.baryU = 1.0f - u - v; // weight for v0
    hit.baryV = u;             // weight for v1
    hit.baryW = v;             // weight for v2

    uint32_t i0 = mesh.indices[triIdx*3+0];
    uint32_t i1 = mesh.indices[triIdx*3+1];
    uint32_t i2 = mesh.indices[triIdx*3+2];

    // Hit position
    hit.position[0] = origin[0] + dir[0] * t;
    hit.position[1] = origin[1] + dir[1] * t;
    hit.position[2] = origin[2] + dir[2] * t;

    // Normal: interpolated if available, otherwise face normal
    if (mesh.hasNormals()) {
        for (int c = 0; c < 3; ++c) {
            hit.normal[c] = hit.baryU * mesh.normals[i0*3+c] +
                            hit.baryV * mesh.normals[i1*3+c] +
                            hit.baryW * mesh.normals[i2*3+c];
        }
        float len = std::sqrt(hit.normal[0]*hit.normal[0] +
                              hit.normal[1]*hit.normal[1] +
                              hit.normal[2]*hit.normal[2]);
        if (len > 1e-8f) {
            hit.normal[0] /= len;
            hit.normal[1] /= len;
            hit.normal[2] /= len;
        }
    } else {
        // Face normal from cross product
        const float* p0 = &mesh.positions[i0*3];
        const float* p1 = &mesh.positions[i1*3];
        const float* p2 = &mesh.positions[i2*3];
        float e1[3] = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
        float e2[3] = {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};
        hit.normal[0] = e1[1]*e2[2] - e1[2]*e2[1];
        hit.normal[1] = e1[2]*e2[0] - e1[0]*e2[2];
        hit.normal[2] = e1[0]*e2[1] - e1[1]*e2[0];
        float len = std::sqrt(hit.normal[0]*hit.normal[0] +
                              hit.normal[1]*hit.normal[1] +
                              hit.normal[2]*hit.normal[2]);
        if (len > 1e-8f) {
            hit.normal[0] /= len;
            hit.normal[1] /= len;
            hit.normal[2] /= len;
        }
    }
}

RayHit raycast(const MeshData& mesh,
               const float* origin, const float* direction,
               float maxDistance) {
    RayHit result;
    if (mesh.empty() || mesh.indices.empty() || !origin || !direction) return result;

    float bestT = (maxDistance > 0.0f) ? maxDistance : std::numeric_limits<float>::max();
    float bestU = 0, bestV = 0;
    uint32_t bestTri = 0;
    bool found = false;

    const size_t triCount = mesh.triangleCount();
    for (size_t tri = 0; tri < triCount; ++tri) {
        uint32_t i0 = mesh.indices[tri*3+0];
        uint32_t i1 = mesh.indices[tri*3+1];
        uint32_t i2 = mesh.indices[tri*3+2];

        float t, u, v;
        if (rayTriangle(origin, direction,
                        &mesh.positions[i0*3], &mesh.positions[i1*3], &mesh.positions[i2*3],
                        t, u, v) && t < bestT) {
            bestT = t;
            bestU = u;
            bestV = v;
            bestTri = (uint32_t)tri;
            found = true;
        }
    }

    if (found) {
        fillHit(result, mesh, origin, direction, bestT, bestU, bestV, bestTri);
    }
    return result;
}

std::vector<RayHit> raycastAll(const MeshData& mesh,
                                const float* origin, const float* direction,
                                float maxDistance) {
    std::vector<RayHit> results;
    if (mesh.empty() || mesh.indices.empty() || !origin || !direction) return results;

    float maxT = (maxDistance > 0.0f) ? maxDistance : std::numeric_limits<float>::max();

    const size_t triCount = mesh.triangleCount();
    for (size_t tri = 0; tri < triCount; ++tri) {
        uint32_t i0 = mesh.indices[tri*3+0];
        uint32_t i1 = mesh.indices[tri*3+1];
        uint32_t i2 = mesh.indices[tri*3+2];

        float t, u, v;
        if (rayTriangle(origin, direction,
                        &mesh.positions[i0*3], &mesh.positions[i1*3], &mesh.positions[i2*3],
                        t, u, v) && t < maxT) {
            RayHit hit;
            fillHit(hit, mesh, origin, direction, t, u, v, (uint32_t)tri);
            results.push_back(hit);
        }
    }

    // Sort by distance
    std::sort(results.begin(), results.end(),
              [](const RayHit& a, const RayHit& b) { return a.distance < b.distance; });

    return results;
}

bool raycastTest(const MeshData& mesh,
                 const float* origin, const float* direction,
                 float maxDistance) {
    if (mesh.empty() || mesh.indices.empty() || !origin || !direction) return false;

    float maxT = (maxDistance > 0.0f) ? maxDistance : std::numeric_limits<float>::max();

    const size_t triCount = mesh.triangleCount();
    for (size_t tri = 0; tri < triCount; ++tri) {
        uint32_t i0 = mesh.indices[tri*3+0];
        uint32_t i1 = mesh.indices[tri*3+1];
        uint32_t i2 = mesh.indices[tri*3+2];

        float t, u, v;
        if (rayTriangle(origin, direction,
                        &mesh.positions[i0*3], &mesh.positions[i1*3], &mesh.positions[i2*3],
                        t, u, v) && t < maxT) {
            return true;
        }
    }
    return false;
}

// Closest point on triangle to a point
static void closestPointOnTriangle(
    const float* p, const float* a, const float* b, const float* c,
    float& outU, float& outV, float& outW, float* outPt) {

    float ab[3] = {b[0]-a[0], b[1]-a[1], b[2]-a[2]};
    float ac[3] = {c[0]-a[0], c[1]-a[1], c[2]-a[2]};
    float ap[3] = {p[0]-a[0], p[1]-a[1], p[2]-a[2]};

    float d1 = ab[0]*ap[0]+ab[1]*ap[1]+ab[2]*ap[2];
    float d2 = ac[0]*ap[0]+ac[1]*ap[1]+ac[2]*ap[2];
    if (d1 <= 0.0f && d2 <= 0.0f) {
        outU = 1; outV = 0; outW = 0;
        outPt[0]=a[0]; outPt[1]=a[1]; outPt[2]=a[2]; return;
    }

    float bp[3] = {p[0]-b[0], p[1]-b[1], p[2]-b[2]};
    float d3 = ab[0]*bp[0]+ab[1]*bp[1]+ab[2]*bp[2];
    float d4 = ac[0]*bp[0]+ac[1]*bp[1]+ac[2]*bp[2];
    if (d3 >= 0.0f && d4 <= d3) {
        outU = 0; outV = 1; outW = 0;
        outPt[0]=b[0]; outPt[1]=b[1]; outPt[2]=b[2]; return;
    }

    float vc = d1*d4 - d3*d2;
    if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
        float v2 = d1 / (d1 - d3);
        outU = 1-v2; outV = v2; outW = 0;
        outPt[0]=a[0]+v2*ab[0]; outPt[1]=a[1]+v2*ab[1]; outPt[2]=a[2]+v2*ab[2]; return;
    }

    float cp[3] = {p[0]-c[0], p[1]-c[1], p[2]-c[2]};
    float d5 = ab[0]*cp[0]+ab[1]*cp[1]+ab[2]*cp[2];
    float d6 = ac[0]*cp[0]+ac[1]*cp[1]+ac[2]*cp[2];
    if (d6 >= 0.0f && d5 <= d6) {
        outU = 0; outV = 0; outW = 1;
        outPt[0]=c[0]; outPt[1]=c[1]; outPt[2]=c[2]; return;
    }

    float vb = d5*d2 - d1*d6;
    if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
        float w2 = d2 / (d2 - d6);
        outU = 1-w2; outV = 0; outW = w2;
        outPt[0]=a[0]+w2*ac[0]; outPt[1]=a[1]+w2*ac[1]; outPt[2]=a[2]+w2*ac[2]; return;
    }

    float va = d3*d6 - d5*d4;
    if (va <= 0.0f && (d4-d3) >= 0.0f && (d5-d6) >= 0.0f) {
        float w2 = (d4-d3) / ((d4-d3)+(d5-d6));
        outU = 0; outV = 1-w2; outW = w2;
        outPt[0]=b[0]+w2*(c[0]-b[0]); outPt[1]=b[1]+w2*(c[1]-b[1]); outPt[2]=b[2]+w2*(c[2]-b[2]); return;
    }

    float denom = 1.0f / (va + vb + vc);
    float v2 = vb * denom;
    float w2 = vc * denom;
    outU = 1-v2-w2; outV = v2; outW = w2;
    outPt[0]=a[0]+ab[0]*v2+ac[0]*w2;
    outPt[1]=a[1]+ab[1]*v2+ac[1]*w2;
    outPt[2]=a[2]+ab[2]*v2+ac[2]*w2;
}

RayHit closestPoint(const MeshData& mesh, const float* point) {
    RayHit result;
    if (mesh.empty() || mesh.indices.empty() || !point) return result;

    float bestDistSq = std::numeric_limits<float>::max();

    const size_t triCount = mesh.triangleCount();
    for (size_t tri = 0; tri < triCount; ++tri) {
        uint32_t i0 = mesh.indices[tri*3+0];
        uint32_t i1 = mesh.indices[tri*3+1];
        uint32_t i2 = mesh.indices[tri*3+2];

        float u, v, w, cp[3];
        closestPointOnTriangle(point,
                               &mesh.positions[i0*3], &mesh.positions[i1*3], &mesh.positions[i2*3],
                               u, v, w, cp);

        float dx = cp[0]-point[0], dy = cp[1]-point[1], dz = cp[2]-point[2];
        float distSq = dx*dx + dy*dy + dz*dz;

        if (distSq < bestDistSq) {
            bestDistSq = distSq;
            result.hit = true;
            result.distance = std::sqrt(distSq);
            result.position[0] = cp[0];
            result.position[1] = cp[1];
            result.position[2] = cp[2];
            result.baryU = u;
            result.baryV = v;
            result.baryW = w;
            result.triangleIndex = (uint32_t)tri;

            // Interpolate normal
            if (mesh.hasNormals()) {
                for (int c2 = 0; c2 < 3; ++c2) {
                    result.normal[c2] = u * mesh.normals[i0*3+c2] +
                                        v * mesh.normals[i1*3+c2] +
                                        w * mesh.normals[i2*3+c2];
                }
                float len = std::sqrt(result.normal[0]*result.normal[0] +
                                      result.normal[1]*result.normal[1] +
                                      result.normal[2]*result.normal[2]);
                if (len > 1e-8f) {
                    result.normal[0] /= len;
                    result.normal[1] /= len;
                    result.normal[2] /= len;
                }
            } else {
                const float* p0 = &mesh.positions[i0*3];
                const float* p1 = &mesh.positions[i1*3];
                const float* p2 = &mesh.positions[i2*3];
                float e1[3] = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
                float e2[3] = {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};
                result.normal[0] = e1[1]*e2[2] - e1[2]*e2[1];
                result.normal[1] = e1[2]*e2[0] - e1[0]*e2[2];
                result.normal[2] = e1[0]*e2[1] - e1[1]*e2[0];
                float len = std::sqrt(result.normal[0]*result.normal[0] +
                                      result.normal[1]*result.normal[1] +
                                      result.normal[2]*result.normal[2]);
                if (len > 1e-8f) {
                    result.normal[0] /= len;
                    result.normal[1] /= len;
                    result.normal[2] /= len;
                }
            }
        }
    }

    return result;
}

} // namespace bromesh
