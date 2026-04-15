#include "bromesh/rigging/auto_rig.h"
#include "bromesh/rigging/skeleton_fit.h"
#include "bromesh/rigging/weight_smooth.h"

#include <cmath>
#include <cstring>

namespace bromesh {

AutoRigResult autoRig(const MeshData& mesh,
                      const RigSpec& spec,
                      const Landmarks& landmarks,
                      const VoxelBindOptions& bindOpts) {
    AutoRigResult r;
    r.missingLandmarks = bromesh::missingLandmarks(spec, landmarks);

    r.skeleton = fitSkeleton(spec, landmarks, mesh);

    // Warn about degenerate (zero-length) bones — usually caused by missing
    // landmarks or bad expressions.
    for (size_t i = 0; i < r.skeleton.bones.size(); ++i) {
        const auto& b = r.skeleton.bones[i];
        float tx = b.localT[0], ty = b.localT[1], tz = b.localT[2];
        float lenSq = tx*tx + ty*ty + tz*tz;
        // Only flag non-root bones with zero offset from parent.
        if (b.parent >= 0 && lenSq < 1e-10f) {
            r.warnings.push_back("zero-length bone: " + b.name);
        }
    }

    r.skin = voxelBindWeights(mesh, r.skeleton, bindOpts);

    // Post-process: Laplacian smoothing + outlier rejection. Safe to call
    // with iterations=0 (just re-normalizes). Driven by the same options
    // struct so callers have a single knob.
    WeightPostProcessOptions smoothOpts;
    smoothOpts.iterations = bindOpts.smoothIterations;
    smoothOpts.alpha      = bindOpts.smoothAlpha;
    smoothOpts.minWeight  = bindOpts.minWeight;
    postProcessWeights(mesh, r.skin, smoothOpts);

    // Sanity check: every vertex should have at least one non-zero weight.
    size_t stride = bindOpts.maxInfluences;
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

} // namespace bromesh
