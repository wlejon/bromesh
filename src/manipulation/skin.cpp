#include "bromesh/manipulation/skin.h"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace bromesh {

void applySkinning(MeshData& mesh, const SkinData& skin,
                   const float* poseMatrices) {
    if (mesh.empty() || skin.boneWeights.empty() || !poseMatrices) return;

    const size_t vCount = mesh.vertexCount();
    if (skin.boneWeights.size() < vCount * 4) return;
    if (skin.boneIndices.size() < vCount * 4) return;

    const bool hasNormals = mesh.hasNormals();

    // Precompute skinning matrices: pose * inverseBindMatrix for each bone
    const size_t boneCount = skin.boneCount;
    std::vector<float> skinMats(boneCount * 16);

    for (size_t b = 0; b < boneCount; ++b) {
        const float* pose = &poseMatrices[b * 16];
        const float* ibm = &skin.inverseBindMatrices[b * 16];
        float* out = &skinMats[b * 16];

        // 4x4 matrix multiply (column-major): out = pose * ibm
        for (int col = 0; col < 4; ++col) {
            for (int row = 0; row < 4; ++row) {
                out[col * 4 + row] =
                    pose[0 * 4 + row] * ibm[col * 4 + 0] +
                    pose[1 * 4 + row] * ibm[col * 4 + 1] +
                    pose[2 * 4 + row] * ibm[col * 4 + 2] +
                    pose[3 * 4 + row] * ibm[col * 4 + 3];
            }
        }
    }

    for (size_t v = 0; v < vCount; ++v) {
        const float* weights = &skin.boneWeights[v * 4];
        const uint32_t* bones = &skin.boneIndices[v * 4];

        float px = mesh.positions[v * 3 + 0];
        float py = mesh.positions[v * 3 + 1];
        float pz = mesh.positions[v * 3 + 2];

        float outPos[3] = {0, 0, 0};
        float outNrm[3] = {0, 0, 0};
        float nx = 0, ny = 0, nz = 0;
        if (hasNormals) {
            nx = mesh.normals[v * 3 + 0];
            ny = mesh.normals[v * 3 + 1];
            nz = mesh.normals[v * 3 + 2];
        }

        for (int j = 0; j < 4; ++j) {
            float w = weights[j];
            if (w <= 0.0f) continue;
            if (bones[j] >= boneCount) continue;

            const float* m = &skinMats[bones[j] * 16];

            // Transform position (affine): m * [px, py, pz, 1]
            outPos[0] += w * (m[0] * px + m[4] * py + m[8]  * pz + m[12]);
            outPos[1] += w * (m[1] * px + m[5] * py + m[9]  * pz + m[13]);
            outPos[2] += w * (m[2] * px + m[6] * py + m[10] * pz + m[14]);

            // Transform normal (linear, no translation): m * [nx, ny, nz, 0]
            if (hasNormals) {
                outNrm[0] += w * (m[0] * nx + m[4] * ny + m[8]  * nz);
                outNrm[1] += w * (m[1] * nx + m[5] * ny + m[9]  * nz);
                outNrm[2] += w * (m[2] * nx + m[6] * ny + m[10] * nz);
            }
        }

        mesh.positions[v * 3 + 0] = outPos[0];
        mesh.positions[v * 3 + 1] = outPos[1];
        mesh.positions[v * 3 + 2] = outPos[2];

        if (hasNormals) {
            float len = std::sqrt(outNrm[0]*outNrm[0] + outNrm[1]*outNrm[1] + outNrm[2]*outNrm[2]);
            if (len > 1e-8f) {
                mesh.normals[v * 3 + 0] = outNrm[0] / len;
                mesh.normals[v * 3 + 1] = outNrm[1] / len;
                mesh.normals[v * 3 + 2] = outNrm[2] / len;
            }
        }
    }
}

void applyMorphTarget(MeshData& mesh, const MorphTarget& morph, float weight) {
    if (mesh.empty() || std::fabs(weight) < 1e-8f) return;

    const size_t vCount = mesh.vertexCount();

    if (morph.deltaPositions.size() >= vCount * 3) {
        for (size_t i = 0; i < vCount * 3; ++i) {
            mesh.positions[i] += morph.deltaPositions[i] * weight;
        }
    }

    if (mesh.hasNormals() && morph.deltaNormals.size() >= vCount * 3) {
        for (size_t i = 0; i < vCount * 3; ++i) {
            mesh.normals[i] += morph.deltaNormals[i] * weight;
        }
        // Re-normalize normals
        for (size_t v = 0; v < vCount; ++v) {
            float* n = &mesh.normals[v * 3];
            float len = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
            if (len > 1e-8f) {
                n[0] /= len; n[1] /= len; n[2] /= len;
            }
        }
    }
}

void normalizeWeights(SkinData& skin) {
    if (skin.boneWeights.empty()) return;

    const size_t vCount = skin.boneWeights.size() / 4;

    for (size_t v = 0; v < vCount; ++v) {
        float* w = &skin.boneWeights[v * 4];
        uint32_t* idx = &skin.boneIndices[v * 4];

        // Zero out near-zero weights
        for (int j = 0; j < 4; ++j) {
            if (w[j] < 1e-6f) {
                w[j] = 0.0f;
                idx[j] = 0;
            }
        }

        // Sort by weight descending (simple bubble for 4 elements)
        for (int a = 0; a < 3; ++a) {
            for (int b = a + 1; b < 4; ++b) {
                if (w[b] > w[a]) {
                    std::swap(w[a], w[b]);
                    std::swap(idx[a], idx[b]);
                }
            }
        }

        // Normalize to sum to 1
        float sum = w[0] + w[1] + w[2] + w[3];
        if (sum > 1e-8f) {
            float inv = 1.0f / sum;
            w[0] *= inv; w[1] *= inv; w[2] *= inv; w[3] *= inv;
        } else {
            // No valid weights — default to bone 0 with full weight
            w[0] = 1.0f; w[1] = 0.0f; w[2] = 0.0f; w[3] = 0.0f;
            idx[0] = 0;
        }
    }
}

} // namespace bromesh
