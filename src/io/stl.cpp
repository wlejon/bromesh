#include "bromesh/io/stl.h"

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>

namespace bromesh {

bool saveSTL(const MeshData& mesh, const std::string& path) {
    if (mesh.empty()) return false;

    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;

    // 80-byte zero header
    char header[80] = {};
    std::fwrite(header, 1, 80, f);

    // Triangle count
    uint32_t triCount = static_cast<uint32_t>(mesh.triangleCount());
    std::fwrite(&triCount, sizeof(uint32_t), 1, f);

    for (uint32_t t = 0; t < triCount; ++t) {
        uint32_t i0 = mesh.indices[t * 3 + 0];
        uint32_t i1 = mesh.indices[t * 3 + 1];
        uint32_t i2 = mesh.indices[t * 3 + 2];

        float v0[3] = { mesh.positions[i0*3+0], mesh.positions[i0*3+1], mesh.positions[i0*3+2] };
        float v1[3] = { mesh.positions[i1*3+0], mesh.positions[i1*3+1], mesh.positions[i1*3+2] };
        float v2[3] = { mesh.positions[i2*3+0], mesh.positions[i2*3+1], mesh.positions[i2*3+2] };

        // Compute face normal: cross(v1-v0, v2-v0)
        float e1[3] = { v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2] };
        float e2[3] = { v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2] };
        float n[3] = {
            e1[1]*e2[2] - e1[2]*e2[1],
            e1[2]*e2[0] - e1[0]*e2[2],
            e1[0]*e2[1] - e1[1]*e2[0]
        };
        float len = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        if (len > 1e-12f) {
            n[0] /= len; n[1] /= len; n[2] /= len;
        }

        std::fwrite(n, sizeof(float), 3, f);
        std::fwrite(v0, sizeof(float), 3, f);
        std::fwrite(v1, sizeof(float), 3, f);
        std::fwrite(v2, sizeof(float), 3, f);

        uint16_t attr = 0;
        std::fwrite(&attr, sizeof(uint16_t), 1, f);
    }

    std::fclose(f);
    return true;
}

MeshData loadSTL(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return {};

    // Skip 80-byte header
    char header[80];
    if (std::fread(header, 1, 80, f) != 80) { std::fclose(f); return {}; }

    // Read triangle count
    uint32_t triCount = 0;
    if (std::fread(&triCount, sizeof(uint32_t), 1, f) != 1) { std::fclose(f); return {}; }

    if (triCount == 0) { std::fclose(f); return {}; }

    MeshData mesh;
    mesh.positions.reserve(triCount * 9);
    mesh.normals.reserve(triCount * 9);
    mesh.indices.reserve(triCount * 3);

    for (uint32_t t = 0; t < triCount; ++t) {
        float normal[3];
        float verts[9]; // 3 vertices * 3 floats
        uint16_t attr;

        if (std::fread(normal, sizeof(float), 3, f) != 3) break;
        if (std::fread(verts, sizeof(float), 9, f) != 9) break;
        if (std::fread(&attr, sizeof(uint16_t), 1, f) != 1) break;

        uint32_t base = static_cast<uint32_t>(mesh.positions.size() / 3);

        for (int v = 0; v < 3; ++v) {
            mesh.positions.push_back(verts[v * 3 + 0]);
            mesh.positions.push_back(verts[v * 3 + 1]);
            mesh.positions.push_back(verts[v * 3 + 2]);
            mesh.normals.push_back(normal[0]);
            mesh.normals.push_back(normal[1]);
            mesh.normals.push_back(normal[2]);
        }

        mesh.indices.push_back(base + 0);
        mesh.indices.push_back(base + 1);
        mesh.indices.push_back(base + 2);
    }

    std::fclose(f);
    return mesh;
}

} // namespace bromesh
