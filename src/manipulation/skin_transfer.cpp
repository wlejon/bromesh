#include "bromesh/manipulation/skin_transfer.h"
#include "bromesh/manipulation/skin.h"
#include "bromesh/analysis/bvh.h"
#include "bromesh/analysis/raycast.h"

#include <algorithm>

namespace bromesh {

SkinData transferSkinWeights(const MeshData& target,
                              const MeshData& source,
                              const SkinData& sourceSkin,
                              float maxDistance) {
    SkinData out;
    out.boneCount = sourceSkin.boneCount;
    out.inverseBindMatrices = sourceSkin.inverseBindMatrices;
    out.boneWeights.assign(target.vertexCount() * 4, 0.0f);
    out.boneIndices.assign(target.vertexCount() * 4, 0u);

    if (source.empty() || target.empty() || sourceSkin.boneWeights.empty()) {
        return out;
    }

    MeshBVH bvh = MeshBVH::build(source);

    // Per-query: collect up to 12 bone contributions (3 source verts x 4
    // slots), merge by bone, keep top 4.
    struct W { uint32_t bone; float weight; };

    for (size_t vi = 0; vi < target.vertexCount(); ++vi) {
        float p[3] = {
            target.positions[vi * 3 + 0],
            target.positions[vi * 3 + 1],
            target.positions[vi * 3 + 2],
        };
        RayHit hit = closestPoint(source, p);
        if (!hit.hit) continue;
        if (maxDistance > 0 && hit.distance > maxDistance) continue;

        uint32_t tri = hit.triangleIndex;
        uint32_t i0 = source.indices[tri * 3 + 0];
        uint32_t i1 = source.indices[tri * 3 + 1];
        uint32_t i2 = source.indices[tri * 3 + 2];
        float b0 = hit.baryU, b1 = hit.baryV, b2 = hit.baryW;

        W bag[12];
        for (int k = 0; k < 12; ++k) { bag[k].bone = 0; bag[k].weight = 0; }

        auto addVert = [&](uint32_t vIdx, float bary, int slot) {
            for (int s = 0; s < 4; ++s) {
                uint32_t b = sourceSkin.boneIndices[vIdx * 4 + s];
                float    w = sourceSkin.boneWeights[vIdx * 4 + s] * bary;
                bag[slot + s].bone = b;
                bag[slot + s].weight = w;
            }
        };
        addVert(i0, b0, 0);
        addVert(i1, b1, 4);
        addVert(i2, b2, 8);

        // Merge duplicate bones
        for (int a = 0; a < 12; ++a) {
            if (bag[a].weight == 0) continue;
            for (int c = a + 1; c < 12; ++c) {
                if (bag[c].weight == 0) continue;
                if (bag[c].bone == bag[a].bone) {
                    bag[a].weight += bag[c].weight;
                    bag[c].weight = 0;
                }
            }
        }

        // Partial sort top 4
        for (int k = 0; k < 4; ++k) {
            int best = k;
            for (int c = k + 1; c < 12; ++c)
                if (bag[c].weight > bag[best].weight) best = c;
            std::swap(bag[k], bag[best]);
        }

        float sum = bag[0].weight + bag[1].weight + bag[2].weight + bag[3].weight;
        if (sum <= 0) continue;
        float inv = 1.0f / sum;
        for (int k = 0; k < 4; ++k) {
            out.boneIndices[vi * 4 + k] = bag[k].bone;
            out.boneWeights[vi * 4 + k] = bag[k].weight * inv;
        }
    }

    return out;
}

} // namespace bromesh
