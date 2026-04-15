#pragma once

#include "bromesh/rigging/rig_spec.h"

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

namespace bromesh {

/// Named 3D points on a mesh. Rigging consumes these as its input (regardless
/// of whether they came from the user clicking in a host app, a heuristic
/// detector, or an ML model). Positions are in the mesh's local space.
struct Landmarks {
    std::unordered_map<std::string, std::array<float,3>> points;

    void set(const std::string& name, float x, float y, float z) {
        points[name] = { x, y, z };
    }

    bool has(const std::string& name) const {
        return points.find(name) != points.end();
    }
};

/// Return the names of landmarks declared by `spec` but missing from `lm`.
/// Empty vector == spec is fully resolvable.
std::vector<std::string> missingLandmarks(const RigSpec& spec,
                                          const Landmarks& lm);

} // namespace bromesh
