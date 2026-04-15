#pragma once

#include <string>
#include <vector>

namespace bromesh {

/// Declarative description of a skeleton template. A RigSpec + a set of
/// landmark positions on a mesh is enough to fit a Skeleton to the mesh and
/// (together with the mesh) compute SkinData.
///
/// A single spec describes any creature: bipeds, quadrupeds, hexapods — the
/// only things that change are the bone list and the landmark schema.
struct RigSpec {
    std::string name;        // "humanoid", "quadruped", ...
    bool symmetric = false;  // enforce left/right mirror during fit

    /// Named mesh landmark the spec expects the caller to provide.
    /// If `mirror` is non-empty, it names the landmark on the opposite side
    /// (e.g. "wrist_L".mirror == "wrist_R"). Used for symmetry enforcement.
    struct LandmarkDecl {
        std::string name;
        std::string mirror;
    };
    std::vector<LandmarkDecl> landmarks;

    /// Bone-position expression, evaluated against the landmark dict to
    /// produce a world-space point. Supported forms:
    ///   "landmark:NAME"              — landmark position
    ///   "mid:A,B"                    — midpoint of landmarks A and B
    ///   "lerp:A,B,t"                 — A + t*(B-A), t in [0,1]
    ///   "offset:A,dx,dy,dz"          — landmark A + constant world offset
    /// The grammar is deliberately tiny. Add forms when real rigs need them;
    /// do not grow it pre-emptively.
    struct BoneDecl {
        std::string name;
        std::string parent;    // empty => root
        std::string head;      // expression for bone head (joint position)
        std::string tail;      // expression for bone tail; empty => auto
        float length = 0.0f;   // fallback tail length when no tail/child
        bool ikChain = false;  // bone is the tip of an IK chain
        bool interior = true;  // project inward off the surface during fit
    };
    std::vector<BoneDecl> bones;

    /// Attachment point on a named bone. Offset is a column-major 4x4.
    struct SocketDecl {
        std::string name;
        std::string bone;
        float offset[16] = {
            1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1
        };
    };
    std::vector<SocketDecl> sockets;
};

/// Parse a RigSpec from JSON text. Returns an empty spec on failure or when
/// JSON support is not compiled in (requires BROMESH_HAS_GLTF for json.hpp).
RigSpec parseRigSpecJSON(const std::string& jsonText);

/// Serialize a RigSpec to JSON text. Returns an empty string when JSON
/// support is not compiled in.
std::string serializeRigSpecJSON(const RigSpec& spec);

/// Load a RigSpec from a JSON file on disk. Returns an empty spec on failure.
RigSpec loadRigSpecFile(const std::string& path);

/// The bundled humanoid spec. Constructed directly in C++ so it works
/// regardless of which optional deps are compiled in.
RigSpec builtinHumanoidSpec();

/// Bundled quadruped spec (dog-like): pelvis-root spine, tail chain, and
/// four bilateral legs. Exists to verify that the RigSpec schema generalizes
/// to non-humanoid topologies without any code changes downstream.
RigSpec builtinQuadrupedSpec();

/// Bundled hexapod spec (insect-like): abdomen-root body and three
/// bilateral leg pairs (front / mid / rear) of three segments each.
RigSpec builtinHexapodSpec();

/// Bundled octopod spec (spider-like): body-root with four bilateral arm
/// pairs of two segments each.
RigSpec builtinOctopodSpec();

/// Look up a bundled spec by name. Accepted values: "humanoid", "quadruped",
/// "hexapod", "octopod". Returns an empty spec for unknown names — callers
/// should check `spec.bones.empty()` before use.
RigSpec builtinRigSpec(const std::string& name);

} // namespace bromesh
