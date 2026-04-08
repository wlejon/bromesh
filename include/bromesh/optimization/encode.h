#pragma once

#include "bromesh/mesh_data.h"
#include <vector>
#include <cstdint>

namespace bromesh {

/// Encoded (compressed) mesh data suitable for GPU streaming.
/// Use decode functions to reconstruct the original data.
struct EncodedMesh {
    std::vector<uint8_t> vertexData;   // encoded vertex buffer
    std::vector<uint8_t> indexData;    // encoded index buffer
    size_t vertexCount;
    size_t vertexSize;                 // bytes per vertex (stride)
    size_t indexCount;
};

/// Encode mesh vertex and index data for efficient streaming/storage.
/// Packs positions (+ optional normals, UVs, colors) into an interleaved
/// vertex buffer and compresses both vertex and index data.
/// Returns empty if meshoptimizer is unavailable.
EncodedMesh encodeMesh(const MeshData& mesh);

/// Decode a previously encoded mesh back to MeshData.
/// Returns empty if decoding fails or meshoptimizer is unavailable.
MeshData decodeMesh(const EncodedMesh& encoded, bool hasNormals, bool hasUVs, bool hasColors);

/// Encode just the index buffer. Returns compressed bytes.
std::vector<uint8_t> encodeIndexBuffer(const std::vector<uint32_t>& indices, size_t vertexCount);

/// Decode index buffer from compressed bytes.
std::vector<uint32_t> decodeIndexBuffer(const std::vector<uint8_t>& data, size_t indexCount);

} // namespace bromesh
