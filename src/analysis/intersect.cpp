#include "bromesh/analysis/intersect.h"

#include <algorithm>
#include <cmath>

namespace bromesh {

// Triangle-triangle intersection test (Möller 1997)
// Tests if two triangles (a0,a1,a2) and (b0,b1,b2) intersect.

static float dot3(const float* a, const float* b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static void cross3(const float* a, const float* b, float* out) {
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
}

static void sub3(const float* a, const float* b, float* out) {
    out[0] = a[0]-b[0]; out[1] = a[1]-b[1]; out[2] = a[2]-b[2];
}

// Compute signed distances of triangle vertices to a plane (normal, d)
static void planeDists(const float* v0, const float* v1, const float* v2,
                       const float* n, float d,
                       float& d0, float& d1, float& d2) {
    d0 = dot3(n, v0) + d;
    d1 = dot3(n, v1) + d;
    d2 = dot3(n, v2) + d;
}

// Compute the interval on the intersection line for a triangle
// Returns false if triangle doesn't cross the line
static bool computeInterval(float d0, float d1, float d2,
                             float p0, float p1, float p2,
                             float& t0, float& t1) {
    // Find the vertex on one side of the plane (the "odd one out")
    // and compute interval endpoints
    if (d0 * d1 > 0.0f) {
        // d2 is the odd one
        t0 = p0 + (p2 - p0) * d0 / (d0 - d2);
        t1 = p1 + (p2 - p1) * d1 / (d1 - d2);
    } else if (d0 * d2 > 0.0f) {
        // d1 is the odd one
        t0 = p0 + (p1 - p0) * d0 / (d0 - d1);
        t1 = p2 + (p1 - p2) * d2 / (d2 - d1);
    } else if (d1 * d2 > 0.0f || d0 != 0.0f) {
        // d0 is the odd one
        t0 = p1 + (p0 - p1) * d1 / (d1 - d0);
        t1 = p2 + (p0 - p2) * d2 / (d2 - d0);
    } else if (d1 != 0.0f) {
        t0 = p0 + (p1 - p0) * d0 / (d0 - d1);
        t1 = p2 + (p1 - p2) * d2 / (d2 - d1);
    } else if (d2 != 0.0f) {
        t0 = p0 + (p2 - p0) * d0 / (d0 - d2);
        t1 = p1 + (p2 - p1) * d1 / (d1 - d2);
    } else {
        // Coplanar — handle conservatively as intersecting
        return false;
    }

    if (t0 > t1) std::swap(t0, t1);
    return true;
}

// ---- Coplanar triangle overlap via 2D separation axis test ----

// Project triangle onto a plane (drop the dominant normal axis) and test overlap
// using the separating axis theorem with edge normals.

static bool coplanarTrianglesOverlap(
    const float* a0, const float* a1, const float* a2,
    const float* b0, const float* b1, const float* b2,
    const float* normal) {

    // Choose projection plane by dropping dominant normal axis
    int dropAxis = 0;
    float absN = std::fabs(normal[0]);
    if (std::fabs(normal[1]) > absN) { dropAxis = 1; absN = std::fabs(normal[1]); }
    if (std::fabs(normal[2]) > absN) { dropAxis = 2; }

    int u = (dropAxis + 1) % 3;
    int v = (dropAxis + 2) % 3;

    float A[3][2] = {{a0[u],a0[v]}, {a1[u],a1[v]}, {a2[u],a2[v]}};
    float B[3][2] = {{b0[u],b0[v]}, {b1[u],b1[v]}, {b2[u],b2[v]}};

    // Separating axis theorem: test all 6 edge normals (3 from A, 3 from B)
    // For each edge, project both triangles onto the perpendicular axis
    for (int pass = 0; pass < 2; ++pass) {
        const float (*tri)[2]  = (pass == 0) ? A : B;
        const float (*other)[2] = (pass == 0) ? B : A;
        for (int i = 0; i < 3; ++i) {
            int j = (i + 1) % 3;
            float ex = tri[j][0] - tri[i][0];
            float ey = tri[j][1] - tri[i][1];
            // Perpendicular axis
            float nx = -ey, ny = ex;

            // Project both triangles onto this axis
            float minA2 = nx*tri[0][0]+ny*tri[0][1], maxA2 = minA2;
            for (int k = 1; k < 3; ++k) {
                float p = nx*tri[k][0]+ny*tri[k][1];
                if (p < minA2) minA2 = p;
                if (p > maxA2) maxA2 = p;
            }
            float minB2 = nx*other[0][0]+ny*other[0][1], maxB2 = minB2;
            for (int k = 1; k < 3; ++k) {
                float p = nx*other[k][0]+ny*other[k][1];
                if (p < minB2) minB2 = p;
                if (p > maxB2) maxB2 = p;
            }

            // Use a small epsilon to avoid false positives on shared edges
            if (maxA2 <= minB2 + 1e-5f || maxB2 <= minA2 + 1e-5f) {
                return false; // Separating axis found
            }
        }
    }

    return true; // No separating axis found => overlap
}

static bool trianglesIntersect(
    const float* a0, const float* a1, const float* a2,
    const float* b0, const float* b1, const float* b2) {

    // Compute plane of triangle A
    float e1[3], e2[3], nA[3];
    sub3(a1, a0, e1);
    sub3(a2, a0, e2);
    cross3(e1, e2, nA);
    float dA = -dot3(nA, a0);

    // Signed distances of B vertices to plane A
    float db0, db1, db2;
    planeDists(b0, b1, b2, nA, dA, db0, db1, db2);

    // Snap near-zero to zero
    if (std::fabs(db0) < 1e-8f) db0 = 0.0f;
    if (std::fabs(db1) < 1e-8f) db1 = 0.0f;
    if (std::fabs(db2) < 1e-8f) db2 = 0.0f;

    // All on same side -> no intersection
    if (db0 > 0 && db1 > 0 && db2 > 0) return false;
    if (db0 < 0 && db1 < 0 && db2 < 0) return false;

    // Compute plane of triangle B
    float f1[3], f2[3], nB[3];
    sub3(b1, b0, f1);
    sub3(b2, b0, f2);
    cross3(f1, f2, nB);
    float dB = -dot3(nB, b0);

    // Signed distances of A vertices to plane B
    float da0, da1, da2;
    planeDists(a0, a1, a2, nB, dB, da0, da1, da2);

    if (std::fabs(da0) < 1e-8f) da0 = 0.0f;
    if (std::fabs(da1) < 1e-8f) da1 = 0.0f;
    if (std::fabs(da2) < 1e-8f) da2 = 0.0f;

    if (da0 > 0 && da1 > 0 && da2 > 0) return false;
    if (da0 < 0 && da1 < 0 && da2 < 0) return false;

    // Compute intersection line direction
    float D[3];
    cross3(nA, nB, D);

    // Project vertices onto intersection line (use largest component of D)
    int maxAxis = 0;
    float maxVal = std::fabs(D[0]);
    if (std::fabs(D[1]) > maxVal) { maxAxis = 1; maxVal = std::fabs(D[1]); }
    if (std::fabs(D[2]) > maxVal) { maxAxis = 2; }

    float pa0 = a0[maxAxis], pa1 = a1[maxAxis], pa2 = a2[maxAxis];
    float pb0 = b0[maxAxis], pb1 = b1[maxAxis], pb2 = b2[maxAxis];

    float tA0, tA1, tB0, tB1;
    bool intA = computeInterval(da0, da1, da2, pa0, pa1, pa2, tA0, tA1);
    bool intB = computeInterval(db0, db1, db2, pb0, pb1, pb2, tB0, tB1);

    if (!intA || !intB) {
        // Coplanar case — use 2D overlap test
        return coplanarTrianglesOverlap(a0, a1, a2, b0, b1, b2, nA);
    }

    // Intervals overlap?
    return tA0 <= tB1 && tB0 <= tA1;
}

// Check if two triangles share any vertex (by index OR by position)
static bool sharesVertex(const uint32_t* a, const uint32_t* b,
                         const float* positions) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (a[i] == b[j]) return true;
            // Also check by position (for meshes with duplicated vertices)
            const float* pa = &positions[a[i]*3];
            const float* pb = &positions[b[j]*3];
            float dx = pa[0]-pb[0], dy = pa[1]-pb[1], dz = pa[2]-pb[2];
            if (dx*dx + dy*dy + dz*dz < 1e-10f) return true;
        }
    }
    return false;
}

// AABB for a triangle
struct TriAABB {
    float min[3], max[3];
};

static TriAABB triBounds(const float* p0, const float* p1, const float* p2) {
    TriAABB b;
    for (int c = 0; c < 3; ++c) {
        b.min[c] = std::min({p0[c], p1[c], p2[c]});
        b.max[c] = std::max({p0[c], p1[c], p2[c]});
    }
    return b;
}

static bool aabbOverlap(const TriAABB& a, const TriAABB& b) {
    for (int c = 0; c < 3; ++c) {
        if (a.min[c] > b.max[c] || b.min[c] > a.max[c]) return false;
    }
    return true;
}

std::vector<TrianglePair> findSelfIntersections(const MeshData& mesh) {
    std::vector<TrianglePair> results;
    if (mesh.empty() || mesh.indices.empty()) return results;

    const size_t triCount = mesh.triangleCount();

    // Precompute AABBs for early rejection
    std::vector<TriAABB> bounds(triCount);
    for (size_t t = 0; t < triCount; ++t) {
        uint32_t i0 = mesh.indices[t*3+0];
        uint32_t i1 = mesh.indices[t*3+1];
        uint32_t i2 = mesh.indices[t*3+2];
        bounds[t] = triBounds(&mesh.positions[i0*3], &mesh.positions[i1*3], &mesh.positions[i2*3]);
    }

    for (size_t a = 0; a < triCount; ++a) {
        for (size_t b = a + 1; b < triCount; ++b) {
            // Skip adjacent triangles
            if (sharesVertex(&mesh.indices[a*3], &mesh.indices[b*3], mesh.positions.data())) continue;

            // AABB early rejection
            if (!aabbOverlap(bounds[a], bounds[b])) continue;

            uint32_t ai0 = mesh.indices[a*3+0], ai1 = mesh.indices[a*3+1], ai2 = mesh.indices[a*3+2];
            uint32_t bi0 = mesh.indices[b*3+0], bi1 = mesh.indices[b*3+1], bi2 = mesh.indices[b*3+2];

            if (trianglesIntersect(
                &mesh.positions[ai0*3], &mesh.positions[ai1*3], &mesh.positions[ai2*3],
                &mesh.positions[bi0*3], &mesh.positions[bi1*3], &mesh.positions[bi2*3])) {
                results.push_back({(uint32_t)a, (uint32_t)b});
            }
        }
    }

    return results;
}

bool hasSelfIntersections(const MeshData& mesh) {
    if (mesh.empty() || mesh.indices.empty()) return false;

    const size_t triCount = mesh.triangleCount();

    std::vector<TriAABB> bounds(triCount);
    for (size_t t = 0; t < triCount; ++t) {
        uint32_t i0 = mesh.indices[t*3+0];
        uint32_t i1 = mesh.indices[t*3+1];
        uint32_t i2 = mesh.indices[t*3+2];
        bounds[t] = triBounds(&mesh.positions[i0*3], &mesh.positions[i1*3], &mesh.positions[i2*3]);
    }

    for (size_t a = 0; a < triCount; ++a) {
        for (size_t b = a + 1; b < triCount; ++b) {
            if (sharesVertex(&mesh.indices[a*3], &mesh.indices[b*3], mesh.positions.data())) continue;
            if (!aabbOverlap(bounds[a], bounds[b])) continue;

            uint32_t ai0 = mesh.indices[a*3+0], ai1 = mesh.indices[a*3+1], ai2 = mesh.indices[a*3+2];
            uint32_t bi0 = mesh.indices[b*3+0], bi1 = mesh.indices[b*3+1], bi2 = mesh.indices[b*3+2];

            if (trianglesIntersect(
                &mesh.positions[ai0*3], &mesh.positions[ai1*3], &mesh.positions[ai2*3],
                &mesh.positions[bi0*3], &mesh.positions[bi1*3], &mesh.positions[bi2*3])) {
                return true;
            }
        }
    }

    return false;
}

bool meshesIntersect(const MeshData& a, const MeshData& b) {
    if (a.empty() || b.empty() || a.indices.empty() || b.indices.empty()) return false;

    const size_t triCountA = a.triangleCount();
    const size_t triCountB = b.triangleCount();

    // Precompute AABBs for mesh B
    std::vector<TriAABB> boundsB(triCountB);
    for (size_t t = 0; t < triCountB; ++t) {
        uint32_t i0 = b.indices[t*3+0], i1 = b.indices[t*3+1], i2 = b.indices[t*3+2];
        boundsB[t] = triBounds(&b.positions[i0*3], &b.positions[i1*3], &b.positions[i2*3]);
    }

    for (size_t ta = 0; ta < triCountA; ++ta) {
        uint32_t ai0 = a.indices[ta*3+0], ai1 = a.indices[ta*3+1], ai2 = a.indices[ta*3+2];
        TriAABB boundsA = triBounds(&a.positions[ai0*3], &a.positions[ai1*3], &a.positions[ai2*3]);

        for (size_t tb = 0; tb < triCountB; ++tb) {
            if (!aabbOverlap(boundsA, boundsB[tb])) continue;

            uint32_t bi0 = b.indices[tb*3+0], bi1 = b.indices[tb*3+1], bi2 = b.indices[tb*3+2];

            if (trianglesIntersect(
                &a.positions[ai0*3], &a.positions[ai1*3], &a.positions[ai2*3],
                &b.positions[bi0*3], &b.positions[bi1*3], &b.positions[bi2*3])) {
                return true;
            }
        }
    }

    return false;
}

} // namespace bromesh
