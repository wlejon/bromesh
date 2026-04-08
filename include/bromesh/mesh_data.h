#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace bromesh {

/// Core mesh output structure. All algorithms produce this.
/// Separate attribute streams for easy TypedArray transfer to JS.
struct MeshData {
    std::vector<float> positions;    // xyz, stride 3
    std::vector<float> normals;      // xyz, stride 3
    std::vector<float> uvs;          // uv,  stride 2 (optional)
    std::vector<float> colors;       // rgba, stride 4 (optional)
    std::vector<uint32_t> indices;

    size_t vertexCount() const { return positions.size() / 3; }
    size_t triangleCount() const { return indices.size() / 3; }
    bool hasNormals() const { return normals.size() == positions.size(); }
    bool hasUVs() const { return uvs.size() / 2 == vertexCount(); }
    bool hasColors() const { return colors.size() / 4 == vertexCount(); }
    bool empty() const { return positions.empty(); }

    void clear() {
        positions.clear();
        normals.clear();
        uvs.clear();
        colors.clear();
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
    float extentX() const { return (max[0] - min[0]) * 0.5f; }
    float extentY() const { return (max[1] - min[1]) * 0.5f; }
    float extentZ() const { return (max[2] - min[2]) * 0.5f; }
};

/// Skinning data for a mesh (optional, populated from glTF).
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

} // namespace bromesh
