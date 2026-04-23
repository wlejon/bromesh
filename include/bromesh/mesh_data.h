#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace bromesh {

/// Core mesh output structure. All algorithms produce this.
/// Separate attribute streams for easy TypedArray transfer to JS.
struct MeshData {
    std::vector<float> positions;    // xyz, stride 3
    std::vector<float> normals;      // xyz, stride 3
    std::vector<float> uvs;          // uv,  stride 2 (optional)
    std::vector<float> colors;       // rgba, stride 4 (optional)
    std::vector<float> tangents;     // xyz + handedness w, stride 4 (optional)
    std::vector<uint32_t> indices;

    size_t vertexCount() const { return positions.size() / 3; }
    size_t triangleCount() const { return indices.size() / 3; }
    bool hasNormals() const { return normals.size() == positions.size(); }
    bool hasUVs() const { return uvs.size() / 2 == vertexCount(); }
    bool hasColors() const { return colors.size() / 4 == vertexCount(); }
    bool hasTangents() const { return tangents.size() / 4 == vertexCount(); }
    bool empty() const { return positions.empty(); }

    void clear() {
        positions.clear();
        normals.clear();
        uvs.clear();
        colors.clear();
        tangents.clear();
        indices.clear();
    }

    /// Reserve space for the given vertex/index counts.
    void reserve(size_t verts, size_t idxCount) {
        positions.reserve(verts * 3);
        normals.reserve(verts * 3);
        indices.reserve(idxCount);
    }
};

/// Axis-aligned bounding box.
struct BBox {
    float min[3] = { 0, 0, 0 };
    float max[3] = { 0, 0, 0 };

    float centerX() const { return (min[0] + max[0]) * 0.5f; }
    float centerY() const { return (min[1] + max[1]) * 0.5f; }
    float centerZ() const { return (min[2] + max[2]) * 0.5f; }

    /// Full size of the box along each axis (max - min). For a [-1, 1] box,
    /// extentX() == 2.0. Matches Three.js / glm naming conventions.
    float extentX() const { return max[0] - min[0]; }
    float extentY() const { return max[1] - min[1]; }
    float extentZ() const { return max[2] - min[2]; }

    /// Half of the full extent (Unity/Unreal "extents" naming). For a [-1, 1]
    /// box, halfExtentX() == 1.0.
    float halfExtentX() const { return (max[0] - min[0]) * 0.5f; }
    float halfExtentY() const { return (max[1] - min[1]) * 0.5f; }
    float halfExtentZ() const { return (max[2] - min[2]) * 0.5f; }
};

/// Skinning data for a mesh (optional, populated from glTF).
/// Per-vertex bone weights and indices; inverse bind matrices are part of
/// Skeleton (kept here too for back-compat with older callers that don't
/// carry a Skeleton alongside).
struct SkinData {
    std::vector<float> boneWeights;    // 4 weights per vertex, stride 4
    std::vector<uint32_t> boneIndices; // 4 indices per vertex, stride 4
    std::vector<float> inverseBindMatrices; // 16 floats (mat4) per bone
    size_t boneCount = 0;
};

/// Morph target: per-vertex deltas.
struct MorphTarget {
    std::string name;
    std::vector<float> deltaPositions; // xyz, stride 3
    std::vector<float> deltaNormals;   // xyz, stride 3 (optional)
};

// --- Skeleton + animation ---------------------------------------------------
//
// A Skeleton is an array of Bones in topological order (parents precede
// children). Bone indices line up with SkinData::boneIndices, so a vertex's
// bone weight index refers to skeleton.bones[i].
//
// Animations are per-bone channels of keyframed T/R/S values. Evaluating an
// Animation at time t produces a pose (local TRS per bone). Composing the
// pose with the skeleton hierarchy + inverse bind matrices gives the joint
// matrices consumed by applySkinning.

struct Bone {
    std::string name;
    int parent = -1;           // index into Skeleton::bones; -1 = root
    float localT[3] = {0,0,0};
    float localR[4] = {0,0,0,1}; // quaternion xyzw
    float localS[3] = {1,1,1};
    float inverseBind[16] = {   // column-major 4x4 mat; identity by default
        1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1
    };
};

/// Named attachment point on the skeleton. The transform chain at runtime is:
/// world(socket) = world(bone) * offset. Use for rigid equipment (weapons,
/// helmets) — no skinning, child follows the bone exactly.
struct Socket {
    std::string name;
    int bone = 0;
    float offset[16] = {
        1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1
    };
};

struct Skeleton {
    std::vector<Bone> bones;
    std::vector<Socket> sockets;

    /// Pre-transform applied to root bones (parent == -1) when computing world
    /// matrices. Column-major 4x4. Identity by default. Used to preserve any
    /// scene-tree ancestry above the skinned joints — e.g. a glTF armature
    /// node with a unit-conversion scale. Keeping it out of the bones' local
    /// TRS means animation keyframes stay in the original (unscaled) space.
    float rootTransform[16] = {
        1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1
    };

    /// Linear bone lookup by name. Returns -1 if not found.
    int findBone(const std::string& boneName) const {
        for (size_t i = 0; i < bones.size(); ++i)
            if (bones[i].name == boneName) return static_cast<int>(i);
        return -1;
    }

    int findSocket(const std::string& socketName) const {
        for (size_t i = 0; i < sockets.size(); ++i)
            if (sockets[i].name == socketName) return static_cast<int>(i);
        return -1;
    }
};

/// Single keyframed channel on a single bone property (T/R/S).
/// Values are packed: translation/scale = 3 floats/key (vec3), rotation = 4
/// floats/key (quaternion xyzw). keyframe count == times.size().
struct AnimChannel {
    enum class Path : uint8_t { Translation, Rotation, Scale };
    enum class Interp : uint8_t { Linear, Step, CubicSpline };

    int boneIndex = -1;
    Path path = Path::Translation;
    Interp interp = Interp::Linear;
    std::vector<float> times;  // seconds, strictly increasing
    std::vector<float> values; // stride = (path == Rotation ? 4 : 3); CubicSpline triples stride
};

struct Animation {
    std::string name;
    float duration = 0.0f; // seconds; max of all channel times
    std::vector<AnimChannel> channels;
};

} // namespace bromesh
