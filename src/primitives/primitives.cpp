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
    int vertCount = (segments + 1) * (rings + 1);
    int idxCount = segments * rings * 6;
    m.reserve(vertCount, idxCount);

    for (int r = 0; r <= rings; r++) {
        float v = (float)r / rings;
        float phi = v * (float)M_PI; // 0 to PI (top to bottom)
        float sinPhi = std::sin(phi);
        float cosPhi = std::cos(phi);

        for (int s = 0; s <= segments; s++) {
            float u = (float)s / segments;
            float theta = u * 2.0f * (float)M_PI;
            float sinTheta = std::sin(theta);
            float cosTheta = std::cos(theta);

            float nx = sinPhi * cosTheta;
            float ny = cosPhi;
            float nz = sinPhi * sinTheta;

            pushVert(m,
                     nx * radius, ny * radius, nz * radius,
                     nx, ny, nz,
                     u, v);
        }
    }

    for (int r = 0; r < rings; r++) {
        for (int s = 0; s < segments; s++) {
            uint32_t a = r * (segments + 1) + s;
            uint32_t b = a + (segments + 1);
            uint32_t c = a + 1;
            uint32_t d = b + 1;

            pushTri(m, a, c, b);
            pushTri(m, c, d, b);
        }
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

    // Top hemisphere: rings/2 bands from pole to equator
    // Cylinder body: 1 band
    // Bottom hemisphere: rings/2 bands from equator to pole
    int halfRings = rings / 2;

    // Total rings of latitude: halfRings (top) + 1 (body top/bottom) + halfRings (bottom)
    // But we build it as a continuous strip

    // We'll build top hemisphere, cylinder sides, then bottom hemisphere
    // Each section shares edge vertices with the next

    // Top hemisphere (from pole down to equator at y = halfHeight)
    for (int r = 0; r <= halfRings; r++) {
        float phi = (float)r / halfRings * ((float)M_PI * 0.5f); // 0 to PI/2
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
                     u, (float)r / (2 * halfRings + 1));
        }
    }

    // Bottom hemisphere (from equator at y = -halfHeight down to pole)
    // Start from r=1 to avoid duplicating the equator ring shared with cylinder
    // Actually, let's add the bottom equator ring and bottom pole
    // The cylinder body connects top hemi equator to bottom hemi equator

    // Bottom equator at y = -halfHeight (same normals as equator = horizontal)
    {
        int r = halfRings + 1;
        for (int s = 0; s <= segments; s++) {
            float u = (float)s / segments;
            float theta = u * 2.0f * (float)M_PI;
            float nx = std::cos(theta);
            float nz = std::sin(theta);
            pushVert(m, nx * radius, -halfHeight, nz * radius, nx, 0, nz,
                     u, (float)r / (2 * halfRings + 1));
        }
    }

    // Bottom hemisphere below -halfHeight
    for (int r = 1; r <= halfRings; r++) {
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
                     u, (float)row / (2 * halfRings + 1));
        }
    }

    // Index all rows as a continuous strip
    int totalRows = 2 * halfRings + 1; // top hemi rows + cylinder + bottom hemi rows
    for (int r = 0; r < totalRows; r++) {
        for (int s = 0; s < segments; s++) {
            uint32_t a = r * (segments + 1) + s;
            uint32_t b = a + (segments + 1);
            uint32_t c = a + 1;
            uint32_t d = b + 1;
            pushTri(m, a, c, b);
            pushTri(m, c, d, b);
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

MeshData heightmapGrid(const float* heights, int gridW, int gridH, float cellSize) {
    if (!heights || gridW < 2 || gridH < 2) return {};

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
            float py = heights[z * gridW + x];
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

    // Compute normals via central differences
    m.normals.resize(vertCount * 3, 0.0f);
    for (int z = 0; z < gridH; z++) {
        for (int x = 0; x < gridW; x++) {
            float hL = heights[z * gridW + (x > 0 ? x - 1 : x)];
            float hR = heights[z * gridW + (x < gridW - 1 ? x + 1 : x)];
            float hD = heights[(z > 0 ? z - 1 : z) * gridW + x];
            float hU = heights[(z < gridH - 1 ? z + 1 : z) * gridW + x];

            float dx = (hR - hL) / ((x > 0 && x < gridW - 1) ? (2.0f * cellSize) : cellSize);
            float dz = (hU - hD) / ((z > 0 && z < gridH - 1) ? (2.0f * cellSize) : cellSize);

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
