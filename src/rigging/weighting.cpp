#include "bromesh/rigging/weighting.h"
#include "bromesh/analysis/bbox.h"

namespace bromesh {

WeightingMethod autoSelectWeightingMethod(const MeshData& mesh) {
    return isManifold(mesh) ? WeightingMethod::BoneHeat
                            : WeightingMethod::VoxelBind;
}

} // namespace bromesh
