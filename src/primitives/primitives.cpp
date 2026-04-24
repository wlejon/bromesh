#include "bromesh/primitives/primitives.h"

#include <cmath>
#include <cstring>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace bromesh {

// Helper: push a vertex (pos + normal + uv)
static void pushVert(MeshData& m, float px, float py, float pz,
                     float nx, float ny, float nz,
                     float u, float v) {
    m.positions.push_back(px);
    m.positions.push_back(py);
    m.positions.push_back(pz);
    m.normals.push_back(nx);
    m.normals.push_back(ny);
    m.normals.push_back(nz);
    m.uvs.push_back(u);
    m.uvs.push_back(v);
}

static void pushTri(MeshData& m, uint32_t a, uint32_t b, uint32_t c) {
    m.indices.push_back(a);
    m.indices.push_back(b);
    m.indices.push_back(c);
}

// ─── box ───────────────────────────────────────────────────────────────────────

MeshData box(float halfW, float halfH, float halfD) {
    MeshData m;
    m.reserve(24, 36);

    // Face data: normal direction, then 4 corner positions
    // Each face has a consistent CCW winding when viewed from outside
    struct Face {
        float nx, ny, nz;
        float v[4][3]; // 4 vertices
        float uv[4][2];
    };

    float w = halfW, h = halfH, d = halfD;

    // +Z face (front)
    auto addFace = [&](float nx, float ny, float nz,
                       float v0x, float v0y, float v0z,
                       float v1x, float v1y, float v1z,
                       float v2x, float v2y, float v2z,
                       float v3x, float v3y, float v3z) {
        uint32_t base = (uint32_t)m.vertexCount();
        pushVert(m, v0x, v0y, v0z, nx, ny, nz, 0, 0);
        pushVert(m, v1x, v1y, v1z, nx, ny, nz, 1, 0);
        pushVert(m, v2x, v2y, v2z, nx, ny, nz, 1, 1);
        pushVert(m, v3x, v3y, v3z, nx, ny, nz, 0, 1);
        pushTri(m, base, base+1, base+2);
        pushTri(m, base, base+2, base+3);
    };

    // +Z (front)
    addFace(0,0,1,  -w,-h, d,  w,-h, d,  w, h, d,  -w, h, d);
    // -Z (back)
    addFace(0,0,-1,  w,-h,-d,  -w,-h,-d,  -w, h,-d,  w, h,-d);
    // +X (right)
    addFace(1,0,0,  w,-h, d,  w,-h,-d,  w, h,-d,  w, h, d);
    // -X (left)
    addFace(-1,0,0,  -w,-h,-d,  -w,-h, d,  -w, h, d,  -w, h,-d);
    // +Y (top)
    addFace(0,1,0,  -w, h, d,  w, h, d,  w, h,-d,  -w, h,-d);
    // -Y (bottom)
    addFace(0,-1,0,  -w,-h,-d,  w,-h,-d,  w,-h, d,  -w,-h, d);

    return m;
}

// ─── sphere ────────────────────────────────────────────────────────────────────

MeshData sphere(float radius, int segments, int rings) {
    MeshData m;
    if (rings < 2 || segments < 3) return m;

    // Pole topology: each pole is a single vertex (UV pinches to u=0.5), not
    // (segments+1) duplicates. The old layout collapsed a whole ring of
    // verts onto each pole position, which broke every topology-aware
    // consumer — half-edge builders saw zero-length edges, CSG saw
    // undefined face normals, OBJ roundtrip dropped the unreferenced
    // duplicates so vertex counts didn't match. A single pole vertex, with
    // a proper fan of triangles under it, is the clean representation.
    const int interiorRings = rings - 1;
    const int vertCount = 2 + interiorRings * (segments + 1);
    const int triCount = 2 * segments + (interiorRings - 1) * 2 * segments;
    m.reserve(vertCount, triCount * 3);

    // Top pole
    pushVert(m, 0.0f, radius, 0.0f, 0.0f, 1.0f, 0.0f, 0.5f, 0.0f);

    // Interior rings r = 1 .. rings-1 (between the poles)
    for (int r = 1; r < rings; r++) {
        float v = (float)r / rings;
        float phi = v * (float)M_PI;
        float sinPhi = std::sin(phi);
        float cosPhi = std::cos(phi);

        for (int s = 0; s <= segments; s++) {
            float u = (float)s / segments;
            float theta = u * 2.0f * (float)M_PI;
            float nx = sinPhi * std::cos(theta);
            float ny = cosPhi;
            float nz = sinPhi * std::sin(theta);

            pushVert(m,
                     nx * radius, ny * radius, nz * radius,
                     nx, ny, nz,
                     u, v);
        }
    }

    // Bottom pole
    pushVert(m, 0.0f, -radius, 0.0f, 0.0f, -1.0f, 0.0f, 0.5f, 1.0f);

    const uint32_t topPole    = 0;
    const uint32_t ring1Start = 1;
    const uint32_t bottomPole = (uint32_t)m.vertexCount() - 1;

    // Top fan
    for (int s = 0; s < segments; s++) {
        uint32_t b = ring1Start + s;
        uint32_t c = ring1Start + s + 1;
        pushTri(m, topPole, c, b);
    }

    // Middle strips (between interior ring r and r+1), r = 1 .. interiorRings-1
    for (int r = 1; r < interiorRings; r++) {
        uint32_t aBase = ring1Start + (r - 1) * (segments + 1);
        uint32_t bBase = ring1Start + r * (segments + 1);
        for (int s = 0; s < segments; s++) {
            uint32_t a = aBase + s;
            uint32_t b = bBase + s;
            uint32_t c = a + 1;
            uint32_t d = b + 1;
            pushTri(m, a, c, b);
            pushTri(m, c, d, b);
        }
    }

    // Bottom fan
    uint32_t lastRingStart = ring1Start + (interiorRings - 1) * (segments + 1);
    for (int s = 0; s < segments; s++) {
        uint32_t a = lastRingStart + s;
        uint32_t c = a + 1;
        pushTri(m, a, c, bottomPole);
    }

    return m;
}

// ─── cylinder ──────────────────────────────────────────────────────────────────

MeshData cylinder(float radius, float halfHeight, int segments) {
    MeshData m;

    // Side vertices: 2 rings of (segments+1) vertices
    // Top cap: 1 center + (segments+1) rim
    // Bottom cap: same
    int sideVerts = (segments + 1) * 2;
    int capVerts = (1 + segments + 1) * 2;
    int sideIdx = segments * 6;
    int capIdx = segments * 3 * 2;
    m.reserve(sideVerts + capVerts, sideIdx + capIdx);

    // Side
    for (int i = 0; i <= 1; i++) {
        float y = (i == 0) ? -halfHeight : halfHeight;
        float v = (float)i;
        for (int s = 0; s <= segments; s++) {
            float u = (float)s / segments;
            float theta = u * 2.0f * (float)M_PI;
            float nx = std::cos(theta);
            float nz = std::sin(theta);
            pushVert(m, nx * radius, y, nz * radius, nx, 0, nz, u, v);
        }
    }

    for (int s = 0; s < segments; s++) {
        uint32_t a = s;
        uint32_t b = s + 1;
        uint32_t c = s + (segments + 1);
        uint32_t d = c + 1;
        pushTri(m, a, c, b);
        pushTri(m, b, c, d);
    }

    // Top cap
    {
        uint32_t centerIdx = (uint32_t)m.vertexCount();
        pushVert(m, 0, halfHeight, 0, 0, 1, 0, 0.5f, 0.5f);
        for (int s = 0; s <= segments; s++) {
            float u = (float)s / segments;
            float theta = u * 2.0f * (float)M_PI;
            float cx = std::cos(theta);
            float cz = std::sin(theta);
            pushVert(m, cx * radius, halfHeight, cz * radius, 0, 1, 0,
                     0.5f + cx * 0.5f, 0.5f + cz * 0.5f);
        }
        for (int s = 0; s < segments; s++) {
            pushTri(m, centerIdx, centerIdx + 2 + s, centerIdx + 1 + s);
        }
    }

    // Bottom cap
    {
        uint32_t centerIdx = (uint32_t)m.vertexCount();
        pushVert(m, 0, -halfHeight, 0, 0, -1, 0, 0.5f, 0.5f);
        for (int s = 0; s <= segments; s++) {
            float u = (float)s / segments;
            float theta = u * 2.0f * (float)M_PI;
            float cx = std::cos(theta);
            float cz = std::sin(theta);
            pushVert(m, cx * radius, -halfHeight, cz * radius, 0, -1, 0,
                     0.5f + cx * 0.5f, 0.5f + cz * 0.5f);
        }
        for (int s = 0; s < segments; s++) {
            pushTri(m, centerIdx, centerIdx + 1 + s, centerIdx + 2 + s);
        }
    }

    return m;
}

// ─── capsule ───────────────────────────────────────────────────────────────────

MeshData capsule(float radius, float halfHeight, int segments, int rings) {
    MeshData m;
    int halfRings = rings / 2;
    if (halfRings < 1 || segments < 3) return m;

    // Pole topology: as in sphere(), each pole is a single vertex — no
    // (segments+1) duplicates that break topology-aware consumers and
    // inflate vertex counts. Interior ring rows are:
    //   top hemi  rings r=1..halfRings     (halfRings rows, each seg+1 verts)
    //   bottom equator                     (1 row, seg+1 verts, at -halfH)
    //   bottom hemi interior r=1..halfRings-1  (halfRings-1 rows)
    // Plus 1 top pole vert and 1 bottom pole vert. 2*halfRings interior
    // ring rows total, plus 2 poles.
    const int interiorRows = 2 * halfRings;
    const int vertCount = 2 + interiorRows * (segments + 1);
    const int triCount = 2 * segments + (interiorRows - 1) * 2 * segments;
    m.reserve(vertCount, triCount * 3);

    const float vDenom = (float)(2 * halfRings + 1);

    // Top pole
    pushVert(m, 0.0f, halfHeight + radius, 0.0f,
             0.0f, 1.0f, 0.0f, 0.5f, 0.0f);

    // Top hemisphere interior rings (r = 1..halfRings). The last is the top
    // equator (phi = pi/2), which the cylinder section connects to below.
    for (int r = 1; r <= halfRings; r++) {
        float phi = (float)r / halfRings * ((float)M_PI * 0.5f);
        float sinPhi = std::sin(phi);
        float cosPhi = std::cos(phi);
        float y = halfHeight + cosPhi * radius;
        for (int s = 0; s <= segments; s++) {
            float u = (float)s / segments;
            float theta = u * 2.0f * (float)M_PI;
            float nx = sinPhi * std::cos(theta);
            float nz = sinPhi * std::sin(theta);
            float ny = cosPhi;
            pushVert(m, nx * radius, y, nz * radius, nx, ny, nz,
                     u, (float)r / vDenom);
        }
    }

    // Bottom equator at y = -halfHeight (row halfRings+1); horizontal normals
    for (int s = 0; s <= segments; s++) {
        float u = (float)s / segments;
        float theta = u * 2.0f * (float)M_PI;
        float nx = std::cos(theta);
        float nz = std::sin(theta);
        pushVert(m, nx * radius, -halfHeight, nz * radius, nx, 0.0f, nz,
                 u, (float)(halfRings + 1) / vDenom);
    }

    // Bottom hemisphere interior rings (r = 1..halfRings-1). The final
    // pole ring of the old layout is replaced by a single pole vertex below.
    for (int r = 1; r < halfRings; r++) {
        float phi = (float)M_PI * 0.5f + (float)r / halfRings * ((float)M_PI * 0.5f);
        float sinPhi = std::sin(phi);
        float cosPhi = std::cos(phi);
        float y = -halfHeight + cosPhi * radius;
        int row = halfRings + 1 + r;
        for (int s = 0; s <= segments; s++) {
            float u = (float)s / segments;
            float theta = u * 2.0f * (float)M_PI;
            float nx = sinPhi * std::cos(theta);
            float nz = sinPhi * std::sin(theta);
            float ny = cosPhi;
            pushVert(m, nx * radius, y, nz * radius, nx, ny, nz,
                     u, (float)row / vDenom);
        }
    }

    // Bottom pole
    pushVert(m, 0.0f, -halfHeight - radius, 0.0f,
             0.0f, -1.0f, 0.0f, 0.5f, 1.0f);

    const uint32_t topPole    = 0;
    const uint32_t ring1Start = 1;
    const uint32_t bottomPole = (uint32_t)m.vertexCount() - 1;

    // Start index of interior row r (r in 1..interiorRows).
    auto rowStart = [&](int r) -> uint32_t {
        return ring1Start + (r - 1) * (segments + 1);
    };

    // Top fan: top pole -> ring 1
    {
        uint32_t ring1 = rowStart(1);
        for (int s = 0; s < segments; s++) {
            uint32_t b = ring1 + s;
            uint32_t c = ring1 + s + 1;
            pushTri(m, topPole, c, b);
        }
    }

    // Interior strips between consecutive interior rows
    for (int r = 1; r < interiorRows; r++) {
        uint32_t aBase = rowStart(r);
        uint32_t bBase = rowStart(r + 1);
        for (int s = 0; s < segments; s++) {
            uint32_t a = aBase + s;
            uint32_t b = bBase + s;
            uint32_t c = a + 1;
            uint32_t d = b + 1;
            pushTri(m, a, c, b);
            pushTri(m, c, d, b);
        }
    }

    // Bottom fan: last interior ring -> bottom pole
    {
        uint32_t lastRing = rowStart(interiorRows);
        for (int s = 0; s < segments; s++) {
            uint32_t a = lastRing + s;
            uint32_t c = a + 1;
            pushTri(m, a, c, bottomPole);
        }
    }

    return m;
}

// ─── plane ─────────────────────────────────────────────────────────────────────

MeshData plane(float halfW, float halfD, int subdivX, int subdivZ) {
    MeshData m;
    int vertsX = subdivX + 1;
    int vertsZ = subdivZ + 1;
    m.reserve(vertsX * vertsZ, subdivX * subdivZ * 6);

    for (int z = 0; z < vertsZ; z++) {
        float v = (float)z / subdivZ;
        float pz = -halfD + v * 2.0f * halfD;
        for (int x = 0; x < vertsX; x++) {
            float u = (float)x / subdivX;
            float px = -halfW + u * 2.0f * halfW;
            pushVert(m, px, 0, pz, 0, 1, 0, u, v);
        }
    }

    for (int z = 0; z < subdivZ; z++) {
        for (int x = 0; x < subdivX; x++) {
            uint32_t a = z * vertsX + x;
            uint32_t b = a + 1;
            uint32_t c = a + vertsX;
            uint32_t d = c + 1;
            pushTri(m, a, c, b);
            pushTri(m, b, c, d);
        }
    }

    return m;
}

// ─── torus ─────────────────────────────────────────────────────────────────────

MeshData torus(float majorRadius, float minorRadius, int majorSegments, int minorSegments) {
    MeshData m;
    int vertCount = (majorSegments + 1) * (minorSegments + 1);
    int idxCount = majorSegments * minorSegments * 6;
    m.reserve(vertCount, idxCount);

    for (int i = 0; i <= majorSegments; i++) {
        float u = (float)i / majorSegments;
        float theta = u * 2.0f * (float)M_PI;
        float cosTheta = std::cos(theta);
        float sinTheta = std::sin(theta);

        for (int j = 0; j <= minorSegments; j++) {
            float v = (float)j / minorSegments;
            float phi = v * 2.0f * (float)M_PI;
            float cosPhi = std::cos(phi);
            float sinPhi = std::sin(phi);

            float px = (majorRadius + minorRadius * cosPhi) * cosTheta;
            float py = minorRadius * sinPhi;
            float pz = (majorRadius + minorRadius * cosPhi) * sinTheta;

            float nx = cosPhi * cosTheta;
            float ny = sinPhi;
            float nz = cosPhi * sinTheta;

            pushVert(m, px, py, pz, nx, ny, nz, u, v);
        }
    }

    for (int i = 0; i < majorSegments; i++) {
        for (int j = 0; j < minorSegments; j++) {
            uint32_t a = i * (minorSegments + 1) + j;
            uint32_t b = a + (minorSegments + 1);
            uint32_t c = a + 1;
            uint32_t d = b + 1;
            pushTri(m, a, c, b);
            pushTri(m, c, d, b);
        }
    }

    return m;
}

// ─── heightmapGrid ─────────────────────────────────────────────────────────────

MeshData heightmapGrid(const float* heights, int gridW, int gridH,
                       float cellSize, int border) {
    if (!heights || gridW < 2 || gridH < 2 || border < 0) return {};

    // Stride through the padded source array.
    int paddedW = gridW + 2 * border;
    auto H = [&](int x, int z) -> float {
        return heights[(z + border) * paddedW + (x + border)];
    };

    MeshData m;
    int vertCount = gridW * gridH;
    int quads = (gridW - 1) * (gridH - 1);
    m.reserve(vertCount, quads * 6);

    float originX = -(gridW - 1) * cellSize * 0.5f;
    float originZ = -(gridH - 1) * cellSize * 0.5f;

    // First pass: positions and UVs
    for (int z = 0; z < gridH; z++) {
        for (int x = 0; x < gridW; x++) {
            float px = originX + x * cellSize;
            float py = H(x, z);
            float pz = originZ + z * cellSize;
            float u = (float)x / (gridW - 1);
            float v = (float)z / (gridH - 1);

            m.positions.push_back(px);
            m.positions.push_back(py);
            m.positions.push_back(pz);
            m.uvs.push_back(u);
            m.uvs.push_back(v);
        }
    }

    // Central-difference normals. With `border >= 1`, the skirt provides the
    // neighbour sample on every side so boundary vertices use the same formula
    // as interior vertices, producing seam-free normals across chunk edges.
    // Without a border we fall back to one-sided differences at the edges.
    m.normals.resize(vertCount * 3, 0.0f);
    for (int z = 0; z < gridH; z++) {
        for (int x = 0; x < gridW; x++) {
            bool hasL = (x > 0) || (border > 0);
            bool hasR = (x < gridW - 1) || (border > 0);
            bool hasD = (z > 0) || (border > 0);
            bool hasU = (z < gridH - 1) || (border > 0);

            float hL = hasL ? H(x - 1, z) : H(x, z);
            float hR = hasR ? H(x + 1, z) : H(x, z);
            float hD = hasD ? H(x, z - 1) : H(x, z);
            float hU = hasU ? H(x, z + 1) : H(x, z);

            float dx = (hR - hL) / ((hasL && hasR) ? (2.0f * cellSize) : cellSize);
            float dz = (hU - hD) / ((hasD && hasU) ? (2.0f * cellSize) : cellSize);

            // Normal = normalize(-dh/dx, 1, -dh/dz)
            float nx = -dx;
            float ny = 1.0f;
            float nz = -dz;
            float len = std::sqrt(nx*nx + ny*ny + nz*nz);
            if (len > 0.0f) { nx /= len; ny /= len; nz /= len; }

            int idx = (z * gridW + x) * 3;
            m.normals[idx]     = nx;
            m.normals[idx + 1] = ny;
            m.normals[idx + 2] = nz;
        }
    }

    // Indices
    for (int z = 0; z < gridH - 1; z++) {
        for (int x = 0; x < gridW - 1; x++) {
            uint32_t a = z * gridW + x;
            uint32_t b = a + 1;
            uint32_t c = a + gridW;
            uint32_t d = c + 1;
            pushTri(m, a, c, b);
            pushTri(m, b, c, d);
        }
    }

    return m;
}

} // namespace bromesh
