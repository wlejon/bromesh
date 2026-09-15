#include "bromesh/manipulation/poly_mesh.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace bromesh {

namespace {

static void cross3(const float a[3], const float b[3], float out[3]) {
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
}

static float dot3(const float a[3], const float b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static float length3(const float a[3]) {
    return std::sqrt(dot3(a, a));
}

static void normalize3(float v[3]) {
    float L = length3(v);
    if (L > 1.0e-20f) { v[0]/=L; v[1]/=L; v[2]/=L; }
    else              { v[0]=v[1]=v[2]=0.0f; }
}

// 2D signed area of triangle a-b-c (positive = CCW).
static float tri2DSignedArea(const float a[2], const float b[2], const float c[2]) {
    return 0.5f * ((b[0]-a[0]) * (c[1]-a[1]) - (b[1]-a[1]) * (c[0]-a[0]));
}

static bool pointInTri2D(const float a[2], const float b[2], const float c[2],
                          const float p[2]) {
    // Barycentric test. Reject on boundary so ear-clipper doesn't pick ears
    // with co-located interior vertices.
    float d1 = tri2DSignedArea(a, b, p);
    float d2 = tri2DSignedArea(b, c, p);
    float d3 = tri2DSignedArea(c, a, p);
    bool hasNeg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    bool hasPos = (d1 > 0) || (d2 > 0) || (d3 > 0);
    return !(hasNeg && hasPos);
}

// Ear-clipping triangulation for a simple CCW polygon in 2D. Writes
// triangles as (i, j, k) index triples referencing uv/idx entries.
// No Steiner points — assumes input polygon is simple and non-self-
// intersecting. Returns false if unable to triangulate (gives up on a
// degenerate/self-intersecting input rather than infinite-looping).
static bool earClip2D(const std::vector<float>& uv,      // stride 2
                      std::vector<uint32_t>& outTris) {
    const int N = (int)(uv.size() / 2);
    if (N < 3) return false;

    // Compute signed area to decide winding.
    float area2 = 0;
    for (int i = 0; i < N; ++i) {
        const float* p = &uv[i*2];
        const float* q = &uv[((i + 1) % N) * 2];
        area2 += p[0]*q[1] - q[0]*p[1];
    }
    const bool ccw = area2 > 0;

    // Working vertex ring as an index list. If CW, walk in reverse so the
    // ear test (requires CCW) works.
    std::vector<int> V(N);
    if (ccw) { for (int i = 0; i < N; ++i) V[i] = i; }
    else     { for (int i = 0; i < N; ++i) V[i] = N - 1 - i; }

    auto at = [&](int idx, float out[2]) {
        out[0] = uv[idx*2 + 0];
        out[1] = uv[idx*2 + 1];
    };

    int remaining = N;
    int guard     = 2 * N;        // degenerate input → bail
    while (remaining > 2) {
        if (--guard < 0) return false;
        bool eared = false;
        for (int i = 0; i < remaining; ++i) {
            int ia = V[(i + remaining - 1) % remaining];
            int ib = V[i];
            int ic = V[(i + 1) % remaining];
            float a[2], b[2], c[2];
            at(ia, a); at(ib, b); at(ic, c);
            // Must be convex (signed area > 0 for CCW).
            if (tri2DSignedArea(a, b, c) <= 0) continue;
            // No other polygon vertex inside.
            bool clean = true;
            for (int j = 0; j < remaining && clean; ++j) {
                int idj = V[j];
                if (idj == ia || idj == ib || idj == ic) continue;
                float p[2]; at(idj, p);
                if (pointInTri2D(a, b, c, p)) clean = false;
            }
            if (!clean) continue;
            outTris.push_back((uint32_t)ia);
            outTris.push_back((uint32_t)ib);
            outTris.push_back((uint32_t)ic);
            V.erase(V.begin() + i);
            --remaining;
            eared = true;
            break;
        }
        if (!eared) return false;   // couldn't find an ear (malformed poly)
    }
    return true;
}

} // namespace

// ===========================================================================
// Tessellation — project each face to its plane, ear-clip in 2D
// ===========================================================================

PolyMesh::Tessellation PolyMesh::tessellate() const {
    Tessellation out;
    if (faces_.empty()) return out;

    // Emit vertices per-face to keep per-face flat normals clean (no sharing
    // across face boundaries). Caller can weld if desired.
    out.positions.reserve(vertices_.size() * 3);
    out.normals.reserve(vertices_.size() * 3);
    out.indices.reserve(faces_.size() * 3);
    out.triToFace.reserve(faces_.size());
    out.triToGroup.reserve(faces_.size());

    for (size_t f = 0; f < faces_.size(); ++f) {
        auto verts = faceVertices((int)f);
        if (verts.size() < 3) continue;

        // Pick plane basis from face normal.
        float n[3]; computeFaceNormal((int)f, n);
        if (dot3(n, n) < 1e-20f) continue;

        // Build u, v on the plane. Choose u ⟂ n avoiding the most-aligned axis.
        float ax = std::abs(n[0]), ay = std::abs(n[1]), az = std::abs(n[2]);
        float tmp[3] = {0,0,0};
        if (ax <= ay && ax <= az)      { tmp[0] = 1; }
        else if (ay <= ax && ay <= az) { tmp[1] = 1; }
        else                           { tmp[2] = 1; }
        float u[3]; cross3(tmp, n, u); normalize3(u);
        float v[3]; cross3(n, u, v);   normalize3(v);

        // Project to 2D.
        std::vector<float> uv(verts.size() * 2);
        for (size_t i = 0; i < verts.size(); ++i) {
            const Vertex& vx = vertices_[verts[i]];
            float p[3] = {vx.x, vx.y, vx.z};
            uv[i*2 + 0] = dot3(p, u);
            uv[i*2 + 1] = dot3(p, v);
        }

        std::vector<uint32_t> tris;
        bool ok = earClip2D(uv, tris);
        if (!ok || tris.empty()) {
            // Fallback: naive fan from vertex 0. Still better than dropping
            // the face entirely — callers can detect degenerate output via
            // validate() and computeFaceNormal.
            tris.clear();
            for (size_t i = 1; i + 1 < verts.size(); ++i) {
                tris.push_back(0);
                tris.push_back((uint32_t)i);
                tris.push_back((uint32_t)(i + 1));
            }
        }

        // Emit positions for this face (per-face dup so shading seams work)
        uint32_t baseIdx = (uint32_t)(out.positions.size() / 3);
        for (int32_t vi : verts) {
            const Vertex& vx = vertices_[vi];
            out.positions.push_back(vx.x);
            out.positions.push_back(vx.y);
            out.positions.push_back(vx.z);
            out.normals.push_back(n[0]);
            out.normals.push_back(n[1]);
            out.normals.push_back(n[2]);
        }
        for (size_t t = 0; t < tris.size(); t += 3) {
            out.indices.push_back(baseIdx + tris[t + 0]);
            out.indices.push_back(baseIdx + tris[t + 1]);
            out.indices.push_back(baseIdx + tris[t + 2]);
            out.triToFace.push_back((int32_t)f);
            out.triToGroup.push_back(faces_[f].group);
        }
    }
    return out;
}

} // namespace bromesh
