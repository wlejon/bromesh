#include "bromesh/uv/projection.h"

#include <cmath>
#include <vector>

namespace bromesh {

static constexpr float PI = 3.14159265358979323846f;

void projectUVs(MeshData& mesh, ProjectionType type, float scale) {
    if (mesh.empty()) return;

    const size_t vCount = mesh.vertexCount();
    mesh.uvs.resize(vCount * 2);

    switch (type) {
        case ProjectionType::PlanarXY: {
            for (size_t i = 0; i < vCount; ++i) {
                mesh.uvs[i * 2 + 0] = mesh.positions[i * 3 + 0] * scale;
                mesh.uvs[i * 2 + 1] = mesh.positions[i * 3 + 1] * scale;
            }
            break;
        }
        case ProjectionType::PlanarXZ: {
            for (size_t i = 0; i < vCount; ++i) {
                mesh.uvs[i * 2 + 0] = mesh.positions[i * 3 + 0] * scale;
                mesh.uvs[i * 2 + 1] = mesh.positions[i * 3 + 2] * scale;
            }
            break;
        }
        case ProjectionType::PlanarYZ: {
            for (size_t i = 0; i < vCount; ++i) {
                mesh.uvs[i * 2 + 0] = mesh.positions[i * 3 + 1] * scale;
                mesh.uvs[i * 2 + 1] = mesh.positions[i * 3 + 2] * scale;
            }
            break;
        }
        case ProjectionType::Cylindrical: {
            for (size_t i = 0; i < vCount; ++i) {
                float x = mesh.positions[i * 3 + 0];
                float y = mesh.positions[i * 3 + 1];
                float z = mesh.positions[i * 3 + 2];
                mesh.uvs[i * 2 + 0] = std::atan2(z, x) / (2.0f * PI);
                mesh.uvs[i * 2 + 1] = y * scale;
            }
            break;
        }
        case ProjectionType::Spherical: {
            for (size_t i = 0; i < vCount; ++i) {
                float x = mesh.positions[i * 3 + 0];
                float y = mesh.positions[i * 3 + 1];
                float z = mesh.positions[i * 3 + 2];
                float r = std::sqrt(x * x + y * y + z * z);
                if (r < 1e-9f) r = 1e-9f;
                mesh.uvs[i * 2 + 0] = std::atan2(z, x) / (2.0f * PI);
                float clamped = y / r;
                if (clamped > 1.0f) clamped = 1.0f;
                if (clamped < -1.0f) clamped = -1.0f;
                mesh.uvs[i * 2 + 1] = std::acos(clamped) / PI;
            }
            break;
        }
        case ProjectionType::Box: {
            // For each vertex, accumulate face normals from adjacent triangles
            // to determine dominant axis, then project onto that plane.
            std::vector<float> avgNormal(vCount * 3, 0.0f);

            const size_t triCount = mesh.triangleCount();
            for (size_t t = 0; t < triCount; ++t) {
                uint32_t i0 = mesh.indices[t * 3 + 0];
                uint32_t i1 = mesh.indices[t * 3 + 1];
                uint32_t i2 = mesh.indices[t * 3 + 2];

                float ax = mesh.positions[i1 * 3 + 0] - mesh.positions[i0 * 3 + 0];
                float ay = mesh.positions[i1 * 3 + 1] - mesh.positions[i0 * 3 + 1];
                float az = mesh.positions[i1 * 3 + 2] - mesh.positions[i0 * 3 + 2];
                float bx = mesh.positions[i2 * 3 + 0] - mesh.positions[i0 * 3 + 0];
                float by = mesh.positions[i2 * 3 + 1] - mesh.positions[i0 * 3 + 1];
                float bz = mesh.positions[i2 * 3 + 2] - mesh.positions[i0 * 3 + 2];

                // Cross product (face normal, not normalized — area-weighted)
                float nx = ay * bz - az * by;
                float ny = az * bx - ax * bz;
                float nz = ax * by - ay * bx;

                for (uint32_t vi : {i0, i1, i2}) {
                    avgNormal[vi * 3 + 0] += nx;
                    avgNormal[vi * 3 + 1] += ny;
                    avgNormal[vi * 3 + 2] += nz;
                }
            }

            // For each vertex, pick dominant axis and project
            for (size_t i = 0; i < vCount; ++i) {
                float anx = std::fabs(avgNormal[i * 3 + 0]);
                float any = std::fabs(avgNormal[i * 3 + 1]);
                float anz = std::fabs(avgNormal[i * 3 + 2]);

                float x = mesh.positions[i * 3 + 0];
                float y = mesh.positions[i * 3 + 1];
                float z = mesh.positions[i * 3 + 2];

                if (anx >= any && anx >= anz) {
                    // Dominant X -> project onto YZ
                    mesh.uvs[i * 2 + 0] = z * scale;
                    mesh.uvs[i * 2 + 1] = y * scale;
                } else if (any >= anx && any >= anz) {
                    // Dominant Y -> project onto XZ
                    mesh.uvs[i * 2 + 0] = x * scale;
                    mesh.uvs[i * 2 + 1] = z * scale;
                } else {
                    // Dominant Z -> project onto XY
                    mesh.uvs[i * 2 + 0] = x * scale;
                    mesh.uvs[i * 2 + 1] = y * scale;
                }
            }
            break;
        }
    }
}

} // namespace bromesh
