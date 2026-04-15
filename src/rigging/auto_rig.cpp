#include "bromesh/rigging/auto_rig.h"
#include "bromesh/rigging/bbw.h"
#include "bromesh/rigging/bone_heat.h"
#include "bromesh/rigging/skeleton_fit.h"
#include "bromesh/rigging/voxel_bind.h"
#include "bromesh/rigging/weight_smooth.h"

#include <cmath>
#include <cstring>

namespace bromesh {

namespace {

SkinData dispatchWeighting(const MeshData& mesh,
                           const Skeleton& skeleton,
                           const WeightingOptions& w,
                           WeightingMethod& chosen) {
    WeightingMethod m = w.method;
    if (m == WeightingMethod::Auto) m = autoSelectWeightingMethod(mesh);
    chosen = m;
    switch (m) {
        case WeightingMethod::BoneHeat: return boneHeatWeights(mesh, skeleton, w.boneHeat);
        case WeightingMethod::BBW:      return bbwWeights(mesh, skeleton, w.bbw);
        case WeightingMethod::VoxelBind:
        default:                        return voxelBindWeights(mesh, skeleton, w.voxel);
    }
}

size_t countInfluences(const WeightingOptions& w, WeightingMethod m) {
    switch (m) {
        case WeightingMethod::BoneHeat: return (size_t)w.boneHeat.maxInfluences;
        case WeightingMethod::BBW:      return (size_t)w.bbw.maxInfluences;
        case WeightingMethod::VoxelBind:
        default:                        return (size_t)w.voxel.maxInfluences;
    }
}

AutoRigResult runAutoRig(const MeshData& mesh,
                         const RigSpec& spec,
                         const Landmarks& landmarks,
                         const WeightingOptions& wopts) {
    AutoRigResult r;
    r.missingLandmarks = bromesh::missingLandmarks(spec, landmarks);
    r.skeleton = fitSkeleton(spec, landmarks, mesh);

    for (size_t i = 0; i < r.skeleton.bones.size(); ++i) {
        const auto& b = r.skeleton.bones[i];
        float tx = b.localT[0], ty = b.localT[1], tz = b.localT[2];
        float lenSq = tx*tx + ty*ty + tz*tz;
        if (b.parent >= 0 && lenSq < 1e-10f) {
            r.warnings.push_back("zero-length bone: " + b.name);
        }
    }

    r.skin = dispatchWeighting(mesh, r.skeleton, wopts, r.methodUsed);

    // Post-process smoothing (safe with 0 iterations; just renormalizes).
    WeightPostProcessOptions smoothOpts;
    smoothOpts.iterations = wopts.smoothIterations;
    smoothOpts.alpha      = wopts.smoothAlpha;
    smoothOpts.minWeight  = wopts.minWeight;
    postProcessWeights(mesh, r.skin, smoothOpts);

    size_t stride = countInfluences(wopts, r.methodUsed);
    size_t orphans = 0;
    for (size_t v = 0; v < mesh.vertexCount(); ++v) {
        bool any = false;
        for (size_t k = 0; k < stride; ++k) {
            if (r.skin.boneWeights[v * stride + k] > 0.0f) { any = true; break; }
        }
        if (!any) ++orphans;
    }
    if (orphans > 0) {
        r.warnings.push_back("orphan vertices without bone influence: "
                             + std::to_string(orphans));
    }
    return r;
}

} // namespace

AutoRigResult autoRig(const MeshData& mesh,
                      const RigSpec& spec,
                      const Landmarks& landmarks,
                      const WeightingOptions& wopts) {
    return runAutoRig(mesh, spec, landmarks, wopts);
}

AutoRigResult autoRig(const MeshData& mesh,
                      const RigSpec& spec,
                      const Landmarks& landmarks,
                      const VoxelBindOptions& bindOpts) {
    WeightingOptions w;
    w.method = WeightingMethod::VoxelBind;
    w.voxel  = bindOpts;
    w.smoothIterations = bindOpts.smoothIterations;
    w.smoothAlpha      = bindOpts.smoothAlpha;
    w.minWeight        = bindOpts.minWeight;
    return runAutoRig(mesh, spec, landmarks, w);
}

} // namespace bromesh
