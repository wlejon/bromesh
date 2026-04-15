#pragma once

#include "bromesh/mesh_data.h"
#include "bromesh/rigging/landmarks.h"
#include "bromesh/rigging/rig_spec.h"
#include "bromesh/rigging/weighting.h"

#include <string>
#include <vector>

namespace bromesh {

struct AutoRigResult {
    Skeleton skeleton;
    SkinData skin;
    std::vector<std::string> missingLandmarks; // required-but-absent landmark names
    std::vector<std::string> warnings;         // soft issues (zero-length bones, orphan verts, ...)
    WeightingMethod methodUsed = WeightingMethod::Auto; // resolved method (never Auto on return)
};

/// One-call driver: validate landmarks, fit skeleton, compute weights.
///
/// Missing landmarks populate `missingLandmarks` in the result but the call
/// still proceeds (affected bones collapse to placeholders). Check the field
/// before consuming the output if strictness matters.
AutoRigResult autoRig(const MeshData& mesh,
                      const RigSpec& spec,
                      const Landmarks& landmarks,
                      const WeightingOptions& wopts = {});

/// Back-compat overload: VoxelBindOptions-only caller. Dispatches to
/// VoxelBind regardless of manifoldness.
AutoRigResult autoRig(const MeshData& mesh,
                      const RigSpec& spec,
                      const Landmarks& landmarks,
                      const VoxelBindOptions& bindOpts);

} // namespace bromesh
