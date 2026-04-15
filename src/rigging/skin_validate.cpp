#include "bromesh/rigging/skin_validate.h"

#include <cmath>

namespace bromesh {

SkinValidation validateSkin(const MeshData& mesh,
                            const SkinData& skin,
                            int influences,
                            float sumTolerance) {
    SkinValidation v;
    v.vertexCount = mesh.vertexCount();
    if (influences <= 0) return v;
    const size_t stride = (size_t)influences;
    if (skin.boneWeights.size() < v.vertexCount * stride) return v;

    for (size_t vi = 0; vi < v.vertexCount; ++vi) {
        float sum = 0.0f;
        int nz = 0;
        bool hasBad = false;
        for (size_t k = 0; k < stride; ++k) {
            float w = skin.boneWeights[vi * stride + k];
            if (std::isnan(w) || w < 0.0f) { hasBad = true; break; }
            sum += w;
            if (w > 0.0f) ++nz;
        }
        if (hasBad) { ++v.nanCount; continue; }
        if (nz == 0) ++v.orphanCount;
        if (nz > v.maxInfluencesObserved) v.maxInfluencesObserved = nz;
        float dev = std::fabs(sum - 1.0f);
        if (dev > v.maxSumDeviation) v.maxSumDeviation = dev;
        if (dev > sumTolerance && nz > 0) ++v.badSumCount;
    }
    return v;
}

} // namespace bromesh
