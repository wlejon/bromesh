#include "bromesh/manipulation/transform.h"
#include "bromesh/isosurface/jit/mir_mesh_kernels.h"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace bromesh {

void transformMesh(MeshData& mesh, const float* m) {
    if (mesh.empty() || !m) return;

    const size_t vCount = mesh.vertexCount();
    const bool hasNormals = mesh.hasNormals();
    const bool hasTangents = mesh.hasTangents();

    const float det3x3 =
        m[0] * (m[5] * m[10] - m[6] * m[9]) -
        m[4] * (m[1] * m[10] - m[2] * m[9]) +
        m[8] * (m[1] * m[6] - m[2] * m[5]);
    const bool flipHandedness = (det3x3 < 0.0f);

    float mats[25];
    std::memcpy(mats, m, 16u * sizeof(float));
    computeNormalMatrix3x3(m, mats + 16);

#if BROMESH_HAS_BRASS_JIT
    if (vCount >= 64 && std::fabs(det3x3) > 1e-12f) {
        static const MirVertexTransformKernel kVertKernel = compileMirVertexTransformKernel();
        if (kVertKernel.cpuValid) {
            const int32_t stride = static_cast<int32_t>(((vCount + 7u) / 8u) * 8u);
            std::vector<float> inSoa(static_cast<size_t>(stride) * 6u, 0.0f);
            std::vector<float> outSoa(static_cast<size_t>(stride) * 6u, 0.0f);

            float* ipx = inSoa.data() + 0 * stride;
            float* ipy = inSoa.data() + 1 * stride;
            float* ipz = inSoa.data() + 2 * stride;
            float* inx = inSoa.data() + 3 * stride;
            float* iny = inSoa.data() + 4 * stride;
            float* inz = inSoa.data() + 5 * stride;
            for (size_t v = 0; v < vCount; ++v) {
                ipx[v] = mesh.positions[v * 3 + 0];
                ipy[v] = mesh.positions[v * 3 + 1];
                ipz[v] = mesh.positions[v * 3 + 2];
                inx[v] = hasNormals ? mesh.normals[v * 3 + 0] : 0.0f;
                iny[v] = hasNormals ? mesh.normals[v * 3 + 1] : 1.0f;
                inz[v] = hasNormals ? mesh.normals[v * 3 + 2] : 0.0f;
            }

            if (kVertKernel.executeCpuSimdSoa(inSoa.data(), outSoa.data(), mats,
                                               static_cast<int32_t>(vCount), stride)) {
                const float* opx = outSoa.data() + 0 * stride;
                const float* opy = outSoa.data() + 1 * stride;
                const float* opz = outSoa.data() + 2 * stride;
                const float* onx = outSoa.data() + 3 * stride;
                const float* ony = outSoa.data() + 4 * stride;
                const float* onz = outSoa.data() + 5 * stride;

                for (size_t v = 0; v < vCount; ++v) {
                    mesh.positions[v * 3 + 0] = opx[v];
                    mesh.positions[v * 3 + 1] = opy[v];
                    mesh.positions[v * 3 + 2] = opz[v];
                    if (hasNormals) {
                        mesh.normals[v * 3 + 0] = onx[v];
                        mesh.normals[v * 3 + 1] = ony[v];
                        mesh.normals[v * 3 + 2] = onz[v];
                    }
                    if (hasTangents) {
                        float tx = mesh.tangents[v * 4 + 0];
                        float ty = mesh.tangents[v * 4 + 1];
                        float tz = mesh.tangents[v * 4 + 2];
                        float tw = mesh.tangents[v * 4 + 3];

                        float ttx = m[0] * tx + m[4] * ty + m[8] * tz;
                        float tty = m[1] * tx + m[5] * ty + m[9] * tz;
                        float ttz = m[2] * tx + m[6] * ty + m[10] * tz;

                        float len = std::sqrt(ttx * ttx + tty * tty + ttz * ttz);
                        if (len > 1e-8f) {
                            mesh.tangents[v * 4 + 0] = ttx / len;
                            mesh.tangents[v * 4 + 1] = tty / len;
                            mesh.tangents[v * 4 + 2] = ttz / len;
                        }
                        mesh.tangents[v * 4 + 3] = flipHandedness ? -tw : tw;
                    }
                }
                return;
            }
        }
    }
#endif

    const float* nmat = mats + 16;
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

            float tnx = nmat[0]*nx + nmat[3]*ny + nmat[6]*nz;
            float tny = nmat[1]*nx + nmat[4]*ny + nmat[7]*nz;
            float tnz = nmat[2]*nx + nmat[5]*ny + nmat[8]*nz;

            float len = std::sqrt(tnx*tnx + tny*tny + tnz*tnz);
            if (len > 1e-8f) {
                mesh.normals[v*3+0] = tnx / len;
                mesh.normals[v*3+1] = tny / len;
                mesh.normals[v*3+2] = tnz / len;
            }
        }

        if (hasTangents) {
            float tx = mesh.tangents[v*4+0];
            float ty = mesh.tangents[v*4+1];
            float tz = mesh.tangents[v*4+2];
            float tw = mesh.tangents[v*4+3];

            float ttx = m[0]*tx + m[4]*ty + m[8]*tz;
            float tty = m[1]*tx + m[5]*ty + m[9]*tz;
            float ttz = m[2]*tx + m[6]*ty + m[10]*tz;

            float len = std::sqrt(ttx*ttx + tty*tty + ttz*ttz);
            if (len > 1e-8f) {
                mesh.tangents[v*4+0] = ttx / len;
                mesh.tangents[v*4+1] = tty / len;
                mesh.tangents[v*4+2] = ttz / len;
            }
            mesh.tangents[v*4+3] = flipHandedness ? -tw : tw;
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
    const bool hasTangents = mesh.hasTangents();

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

        if (hasTangents) {
            float tx = mesh.tangents[v*4+0] * sx;
            float ty = mesh.tangents[v*4+1] * sy;
            float tz = mesh.tangents[v*4+2] * sz;
            float len = std::sqrt(tx*tx + ty*ty + tz*tz);
            if (len > 1e-8f) {
                mesh.tangents[v*4+0] = tx / len;
                mesh.tangents[v*4+1] = ty / len;
                mesh.tangents[v*4+2] = tz / len;
            }
        }
    }

    // If any axis is negative, winding flips for an odd number of negative axes
    int negCount = (sx < 0 ? 1 : 0) + (sy < 0 ? 1 : 0) + (sz < 0 ? 1 : 0);
    if (negCount % 2 == 1) {
        for (size_t t = 0; t < mesh.triangleCount(); ++t) {
            std::swap(mesh.indices[t*3+1], mesh.indices[t*3+2]);
        }
        if (hasTangents) {
            for (size_t v = 0; v < vCount; ++v) {
                mesh.tangents[v*4+3] = -mesh.tangents[v*4+3];
            }
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
    const bool hasTangents = mesh.hasTangents();
    for (size_t v = 0; v < vCount; ++v) {
        mesh.positions[v*3 + axis] = -mesh.positions[v*3 + axis];
        if (mesh.hasNormals()) {
            mesh.normals[v*3 + axis] = -mesh.normals[v*3 + axis];
        }
        if (hasTangents) {
            mesh.tangents[v*4 + axis] = -mesh.tangents[v*4 + axis];
            mesh.tangents[v*4 + 3] = -mesh.tangents[v*4 + 3];
        }
    }

    // Reverse winding order to preserve face orientation
    for (size_t t = 0; t < mesh.triangleCount(); ++t) {
        std::swap(mesh.indices[t*3+1], mesh.indices[t*3+2]);
    }
}

bromath::Vec3 centerMesh(MeshData& mesh) {
    if (mesh.empty()) return {0.0f, 0.0f, 0.0f};

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

    translateMesh(mesh, -cx, -cy, -cz);
    return {cx, cy, cz};
}

} // namespace bromesh
