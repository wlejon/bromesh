#include "bromesh/manipulation/transform.h"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace bromesh {

void transformMesh(MeshData& mesh, const float* m) {
    if (mesh.empty() || !m) return;

    const size_t vCount = mesh.vertexCount();
    const bool hasNormals = mesh.hasNormals();

    // Compute inverse-transpose of upper 3x3 for normal transformation.
    // For a general 3x3 A, inv(A)^T = adj(A)^T / det(A).
    // We only need the direction so we can skip dividing by det (and re-normalize).
    float n00 = m[5]*m[10] - m[6]*m[9];
    float n01 = m[6]*m[8]  - m[4]*m[10];
    float n02 = m[4]*m[9]  - m[5]*m[8];
    float n10 = m[2]*m[9]  - m[1]*m[10];
    float n11 = m[0]*m[10] - m[2]*m[8];
    float n12 = m[1]*m[8]  - m[0]*m[9];
    float n20 = m[1]*m[6]  - m[2]*m[5];
    float n21 = m[2]*m[4]  - m[0]*m[6];
    float n22 = m[0]*m[5]  - m[1]*m[4];

    for (size_t v = 0; v < vCount; ++v) {
        float x = mesh.positions[v*3+0];
        float y = mesh.positions[v*3+1];
        float z = mesh.positions[v*3+2];

        // Column-major: pos = M * [x,y,z,1]
        mesh.positions[v*3+0] = m[0]*x + m[4]*y + m[8]*z  + m[12];
        mesh.positions[v*3+1] = m[1]*x + m[5]*y + m[9]*z  + m[13];
        mesh.positions[v*3+2] = m[2]*x + m[6]*y + m[10]*z + m[14];

        if (hasNormals) {
            float nx = mesh.normals[v*3+0];
            float ny = mesh.normals[v*3+1];
            float nz = mesh.normals[v*3+2];

            float tnx = n00*nx + n10*ny + n20*nz;
            float tny = n01*nx + n11*ny + n21*nz;
            float tnz = n02*nx + n12*ny + n22*nz;

            float len = std::sqrt(tnx*tnx + tny*tny + tnz*tnz);
            if (len > 1e-8f) {
                mesh.normals[v*3+0] = tnx / len;
                mesh.normals[v*3+1] = tny / len;
                mesh.normals[v*3+2] = tnz / len;
            }
        }
    }
}

void translateMesh(MeshData& mesh, float dx, float dy, float dz) {
    const size_t vCount = mesh.vertexCount();
    for (size_t v = 0; v < vCount; ++v) {
        mesh.positions[v*3+0] += dx;
        mesh.positions[v*3+1] += dy;
        mesh.positions[v*3+2] += dz;
    }
}

void scaleMesh(MeshData& mesh, float sx, float sy, float sz) {
    const size_t vCount = mesh.vertexCount();
    const bool hasNormals = mesh.hasNormals();

    for (size_t v = 0; v < vCount; ++v) {
        mesh.positions[v*3+0] *= sx;
        mesh.positions[v*3+1] *= sy;
        mesh.positions[v*3+2] *= sz;

        if (hasNormals) {
            // Scale normals by inverse scale and re-normalize
            float nx = mesh.normals[v*3+0] / sx;
            float ny = mesh.normals[v*3+1] / sy;
            float nz = mesh.normals[v*3+2] / sz;
            float len = std::sqrt(nx*nx + ny*ny + nz*nz);
            if (len > 1e-8f) {
                mesh.normals[v*3+0] = nx / len;
                mesh.normals[v*3+1] = ny / len;
                mesh.normals[v*3+2] = nz / len;
            }
        }
    }

    // If any axis is negative, winding flips for an odd number of negative axes
    int negCount = (sx < 0 ? 1 : 0) + (sy < 0 ? 1 : 0) + (sz < 0 ? 1 : 0);
    if (negCount % 2 == 1) {
        for (size_t t = 0; t < mesh.triangleCount(); ++t) {
            std::swap(mesh.indices[t*3+1], mesh.indices[t*3+2]);
        }
    }
}

void scaleMesh(MeshData& mesh, float s) {
    scaleMesh(mesh, s, s, s);
}

void rotateMesh(MeshData& mesh, float ax, float ay, float az, float angleRadians) {
    // Normalize axis
    float len = std::sqrt(ax*ax + ay*ay + az*az);
    if (len < 1e-8f) return;
    ax /= len; ay /= len; az /= len;

    // Build rotation matrix (Rodrigues' formula), column-major
    float c = std::cos(angleRadians);
    float s = std::sin(angleRadians);
    float t = 1.0f - c;

    float m[16] = {
        t*ax*ax + c,      t*ax*ay + s*az,   t*ax*az - s*ay,   0,
        t*ax*ay - s*az,   t*ay*ay + c,      t*ay*az + s*ax,   0,
        t*ax*az + s*ay,   t*ay*az - s*ax,   t*az*az + c,      0,
        0,                0,                0,                1
    };

    transformMesh(mesh, m);
}

void mirrorMesh(MeshData& mesh, int axis) {
    if (axis < 0 || axis > 2) return;

    const size_t vCount = mesh.vertexCount();
    for (size_t v = 0; v < vCount; ++v) {
        mesh.positions[v*3 + axis] = -mesh.positions[v*3 + axis];
        if (mesh.hasNormals()) {
            mesh.normals[v*3 + axis] = -mesh.normals[v*3 + axis];
        }
    }

    // Reverse winding order to preserve face orientation
    for (size_t t = 0; t < mesh.triangleCount(); ++t) {
        std::swap(mesh.indices[t*3+1], mesh.indices[t*3+2]);
    }
}

void centerMesh(MeshData& mesh, float* outCenter) {
    if (mesh.empty()) return;

    const size_t vCount = mesh.vertexCount();
    float bmin[3] = {mesh.positions[0], mesh.positions[1], mesh.positions[2]};
    float bmax[3] = {mesh.positions[0], mesh.positions[1], mesh.positions[2]};
    for (size_t v = 1; v < vCount; ++v) {
        for (int c = 0; c < 3; ++c) {
            bmin[c] = std::min(bmin[c], mesh.positions[v*3+c]);
            bmax[c] = std::max(bmax[c], mesh.positions[v*3+c]);
        }
    }

    float cx = (bmin[0] + bmax[0]) * 0.5f;
    float cy = (bmin[1] + bmax[1]) * 0.5f;
    float cz = (bmin[2] + bmax[2]) * 0.5f;

    if (outCenter) {
        outCenter[0] = cx;
        outCenter[1] = cy;
        outCenter[2] = cz;
    }

    translateMesh(mesh, -cx, -cy, -cz);
}

} // namespace bromesh
