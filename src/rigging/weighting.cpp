#include "bromesh/rigging/weighting.h"
#include "bromesh/analysis/bbox.h"

#include <cstring>

namespace bromesh {

WeightingMethod autoSelectWeightingMethod(const MeshData& mesh) {
    return isManifold(mesh) ? WeightingMethod::BoneHeat
                            : WeightingMethod::VoxelBind;
}

const char* weightingMethodName(WeightingMethod m) {
    switch (m) {
        case WeightingMethod::Auto:      return "auto";
        case WeightingMethod::VoxelBind: return "voxel";
        case WeightingMethod::BoneHeat:  return "boneHeat";
        case WeightingMethod::BBW:       return "bbw";
    }
    return "auto";
}

WeightingMethod parseWeightingMethod(const char* name) {
    if (!name) return WeightingMethod::Auto;
    if (std::strcmp(name, "voxel")    == 0) return WeightingMethod::VoxelBind;
    if (std::strcmp(name, "boneHeat") == 0) return WeightingMethod::BoneHeat;
    if (std::strcmp(name, "bbw")      == 0) return WeightingMethod::BBW;
    return WeightingMethod::Auto;
}

} // namespace bromesh
