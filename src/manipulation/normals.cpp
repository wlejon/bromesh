#include "bromesh/manipulation/normals.h"
#include <cmath>
#include <cstring>

namespace bromesh {

void computeNormals(MeshData& mesh) {
    if (mesh.positions.empty() || mesh.indices.empty()) return;

    size_t vertCount = mesh.vertexCount();
    mesh.normals.assign(vertCount * 3, 0.0f);

    // Accumulate area-weighted face normals per vertex
    size_t triCount = mesh.triangleCount();
    for (size_t t = 0; t < triCount; ++t) {
        uint32_t i0 = mesh.indices[t * 3 + 0];
        uint32_t i1 = mesh.indices[t * 3 + 1];
        uint32_t i2 = mesh.indices[t * 3 + 2];

        float ax = mesh.positions[i0 * 3 + 0];
        float ay = mesh.positions[i0 * 3 + 1];
        float az = mesh.positions[i0 * 3 + 2];
        float bx = mesh.positions[i1 * 3 + 0];
        float by = mesh.positions[i1 * 3 + 1];
        float bz = mesh.positions[i1 * 3 + 2];
        float cx = mesh.positions[i2 * 3 + 0];
        float cy = mesh.positions[i2 * 3 + 1];
        float cz = mesh.positions[i2 * 3 + 2];

        // Edge vectors
        float e1x = bx - ax, e1y = by - ay, e1z = bz - az;
        float e2x = cx - ax, e2y = cy - ay, e2z = cz - az;

        // Cross product (magnitude = 2 * triangle area)
        float nx = e1y * e2z - e1z * e2y;
        float ny = e1z * e2x - e1x * e2z;
        float nz = e1x * e2y - e1y * e2x;

        // Accumulate (area-weighted since we don't normalize the face normal)
        for (uint32_t idx : {i0, i1, i2}) {
            mesh.normals[idx * 3 + 0] += nx;
            mesh.normals[idx * 3 + 1] += ny;
            mesh.normals[idx * 3 + 2] += nz;
        }
    }

    // Normalize
    for (size_t v = 0; v < vertCount; ++v) {
        float nx = mesh.normals[v * 3 + 0];
        float ny = mesh.normals[v * 3 + 1];
        float nz = mesh.normals[v * 3 + 2];
        float len = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (len > 1e-8f) {
            float inv = 1.0f / len;
            mesh.normals[v * 3 + 0] *= inv;
            mesh.normals[v * 3 + 1] *= inv;
            mesh.normals[v * 3 + 2] *= inv;
        }
    }
}

MeshData computeFlatNormals(const MeshData& mesh) {
    if (mesh.positions.empty() || mesh.indices.empty()) return {};

    MeshData out;
    size_t triCount = mesh.triangleCount();
    size_t newVertCount = triCount * 3;

    out.positions.resize(newVertCount * 3);
    out.normals.resize(newVertCount * 3);
    out.indices.resize(newVertCount);

    bool hasUVs = mesh.hasUVs();
    if (hasUVs) out.uvs.resize(newVertCount * 2);

    bool hasColors = mesh.hasColors();
    if (hasColors) out.colors.resize(newVertCount * 4);

    for (size_t t = 0; t < triCount; ++t) {
        uint32_t i0 = mesh.indices[t * 3 + 0];
        uint32_t i1 = mesh.indices[t * 3 + 1];
        uint32_t i2 = mesh.indices[t * 3 + 2];

        // Copy positions
        for (int k = 0; k < 3; ++k) {
            out.positions[(t * 3 + 0) * 3 + k] = mesh.positions[i0 * 3 + k];
            out.positions[(t * 3 + 1) * 3 + k] = mesh.positions[i1 * 3 + k];
            out.positions[(t * 3 + 2) * 3 + k] = mesh.positions[i2 * 3 + k];
        }

        // Compute face normal
        float ax = mesh.positions[i0 * 3], ay = mesh.positions[i0 * 3 + 1], az = mesh.positions[i0 * 3 + 2];
        float bx = mesh.positions[i1 * 3], by = mesh.positions[i1 * 3 + 1], bz = mesh.positions[i1 * 3 + 2];
        float cx = mesh.positions[i2 * 3], cy = mesh.positions[i2 * 3 + 1], cz = mesh.positions[i2 * 3 + 2];

        float e1x = bx - ax, e1y = by - ay, e1z = bz - az;
        float e2x = cx - ax, e2y = cy - ay, e2z = cz - az;
        float nx = e1y * e2z - e1z * e2y;
        float ny = e1z * e2x - e1x * e2z;
        float nz = e1x * e2y - e1y * e2x;
        float len = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (len > 1e-8f) {
            float inv = 1.0f / len;
            nx *= inv; ny *= inv; nz *= inv;
        }

        for (int v = 0; v < 3; ++v) {
            out.normals[(t * 3 + v) * 3 + 0] = nx;
            out.normals[(t * 3 + v) * 3 + 1] = ny;
            out.normals[(t * 3 + v) * 3 + 2] = nz;
        }

        // Copy UVs
        if (hasUVs) {
            uint32_t ids[3] = {i0, i1, i2};
            for (int v = 0; v < 3; ++v) {
                out.uvs[(t * 3 + v) * 2 + 0] = mesh.uvs[ids[v] * 2 + 0];
                out.uvs[(t * 3 + v) * 2 + 1] = mesh.uvs[ids[v] * 2 + 1];
            }
        }

        // Copy colors
        if (hasColors) {
            uint32_t ids[3] = {i0, i1, i2};
            for (int v = 0; v < 3; ++v) {
                for (int c = 0; c < 4; ++c)
                    out.colors[(t * 3 + v) * 4 + c] = mesh.colors[ids[v] * 4 + c];
            }
        }

        // Sequential indices
        out.indices[t * 3 + 0] = static_cast<uint32_t>(t * 3 + 0);
        out.indices[t * 3 + 1] = static_cast<uint32_t>(t * 3 + 1);
        out.indices[t * 3 + 2] = static_cast<uint32_t>(t * 3 + 2);
    }

    return out;
}

std::vector<float> computeTangents(const MeshData& mesh) {
    if (!mesh.hasUVs() || mesh.positions.empty() || mesh.indices.empty())
        return {};

    size_t vertCount = mesh.vertexCount();
    std::vector<float> tan1(vertCount * 3, 0.0f); // tangent accumulator
    std::vector<float> tan2(vertCount * 3, 0.0f); // bitangent accumulator

    size_t triCount = mesh.triangleCount();
    for (size_t t = 0; t < triCount; ++t) {
        uint32_t i0 = mesh.indices[t * 3 + 0];
        uint32_t i1 = mesh.indices[t * 3 + 1];
        uint32_t i2 = mesh.indices[t * 3 + 2];

        float p0x = mesh.positions[i0 * 3], p0y = mesh.positions[i0 * 3 + 1], p0z = mesh.positions[i0 * 3 + 2];
        float p1x = mesh.positions[i1 * 3], p1y = mesh.positions[i1 * 3 + 1], p1z = mesh.positions[i1 * 3 + 2];
        float p2x = mesh.positions[i2 * 3], p2y = mesh.positions[i2 * 3 + 1], p2z = mesh.positions[i2 * 3 + 2];

        float u0 = mesh.uvs[i0 * 2], v0 = mesh.uvs[i0 * 2 + 1];
        float u1 = mesh.uvs[i1 * 2], v1 = mesh.uvs[i1 * 2 + 1];
        float u2 = mesh.uvs[i2 * 2], v2 = mesh.uvs[i2 * 2 + 1];

        float dx1 = p1x - p0x, dy1 = p1y - p0y, dz1 = p1z - p0z;
        float dx2 = p2x - p0x, dy2 = p2y - p0y, dz2 = p2z - p0z;
        float du1 = u1 - u0, dv1 = v1 - v0;
        float du2 = u2 - u0, dv2 = v2 - v0;

        float det = du1 * dv2 - du2 * dv1;
        if (std::fabs(det) < 1e-8f) continue;
        float r = 1.0f / det;

        float tx = (dv2 * dx1 - dv1 * dx2) * r;
        float ty = (dv2 * dy1 - dv1 * dy2) * r;
        float tz = (dv2 * dz1 - dv1 * dz2) * r;

        float bx = (du1 * dx2 - du2 * dx1) * r;
        float by = (du1 * dy2 - du2 * dy1) * r;
        float bz = (du1 * dz2 - du2 * dz1) * r;

        for (uint32_t idx : {i0, i1, i2}) {
            tan1[idx * 3 + 0] += tx;
            tan1[idx * 3 + 1] += ty;
            tan1[idx * 3 + 2] += tz;
            tan2[idx * 3 + 0] += bx;
            tan2[idx * 3 + 1] += by;
            tan2[idx * 3 + 2] += bz;
        }
    }

    // Gram-Schmidt orthogonalize and compute handedness
    std::vector<float> tangents(vertCount * 4);
    for (size_t v = 0; v < vertCount; ++v) {
        float nx = mesh.normals[v * 3 + 0];
        float ny = mesh.normals[v * 3 + 1];
        float nz = mesh.normals[v * 3 + 2];
        float tx = tan1[v * 3 + 0];
        float ty = tan1[v * 3 + 1];
        float tz = tan1[v * 3 + 2];

        // Gram-Schmidt: t' = t - n * dot(n, t)
        float dot = nx * tx + ny * ty + nz * tz;
        float ox = tx - nx * dot;
        float oy = ty - ny * dot;
        float oz = tz - nz * dot;
        float len = std::sqrt(ox * ox + oy * oy + oz * oz);
        if (len > 1e-8f) {
            float inv = 1.0f / len;
            ox *= inv; oy *= inv; oz *= inv;
        }

        // Handedness: sign of dot(cross(n, t), tan2)
        float cx = ny * tz - nz * ty;
        float cy = nz * tx - nx * tz;
        float cz = nx * ty - ny * tx;
        float w = (cx * tan2[v * 3] + cy * tan2[v * 3 + 1] + cz * tan2[v * 3 + 2]) < 0.0f ? -1.0f : 1.0f;

        tangents[v * 4 + 0] = ox;
        tangents[v * 4 + 1] = oy;
        tangents[v * 4 + 2] = oz;
        tangents[v * 4 + 3] = w;
    }

    return tangents;
}

} // namespace bromesh
