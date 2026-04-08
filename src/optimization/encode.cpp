#include "bromesh/optimization/encode.h"

#ifdef BROMESH_HAS_MESHOPTIMIZER
#include <meshoptimizer.h>
#include <cstring>
#endif

namespace bromesh {

EncodedMesh encodeMesh(const MeshData& mesh) {
#ifdef BROMESH_HAS_MESHOPTIMIZER
    if (mesh.empty() || mesh.indices.empty()) return {};

    size_t vertexCount = mesh.vertexCount();
    size_t indexCount = mesh.indices.size();

    // Compute vertex stride: positions(12) + normals(12) + uvs(8) + colors(16)
    size_t stride = 12; // positions always present
    if (mesh.hasNormals()) stride += 12;
    if (mesh.hasUVs()) stride += 8;
    if (mesh.hasColors()) stride += 16;

    // Interleave vertex data
    std::vector<uint8_t> interleavedVertices(vertexCount * stride);
    for (size_t v = 0; v < vertexCount; ++v) {
        uint8_t* dst = interleavedVertices.data() + v * stride;
        size_t off = 0;

        std::memcpy(dst + off, &mesh.positions[v * 3], 12);
        off += 12;

        if (mesh.hasNormals()) {
            std::memcpy(dst + off, &mesh.normals[v * 3], 12);
            off += 12;
        }
        if (mesh.hasUVs()) {
            std::memcpy(dst + off, &mesh.uvs[v * 2], 8);
            off += 8;
        }
        if (mesh.hasColors()) {
            std::memcpy(dst + off, &mesh.colors[v * 4], 16);
            off += 16;
        }
    }

    // Encode vertex buffer
    size_t vbBound = meshopt_encodeVertexBufferBound(vertexCount, stride);
    std::vector<uint8_t> encodedVB(vbBound);
    size_t vbSize = meshopt_encodeVertexBuffer(
        encodedVB.data(), encodedVB.size(),
        interleavedVertices.data(), vertexCount, stride
    );
    encodedVB.resize(vbSize);

    // Encode index buffer
    size_t ibBound = meshopt_encodeIndexBufferBound(indexCount, vertexCount);
    std::vector<uint8_t> encodedIB(ibBound);
    size_t ibSize = meshopt_encodeIndexBuffer(
        encodedIB.data(), encodedIB.size(),
        mesh.indices.data(), indexCount
    );
    encodedIB.resize(ibSize);

    EncodedMesh result;
    result.vertexData = std::move(encodedVB);
    result.indexData = std::move(encodedIB);
    result.vertexCount = vertexCount;
    result.vertexSize = stride;
    result.indexCount = indexCount;
    return result;
#else
    (void)mesh;
    return {};
#endif
}

MeshData decodeMesh(const EncodedMesh& encoded, bool hasNormals, bool hasUVs, bool hasColors) {
#ifdef BROMESH_HAS_MESHOPTIMIZER
    if (encoded.vertexData.empty() || encoded.indexData.empty()) return {};

    size_t vertexCount = encoded.vertexCount;
    size_t stride = encoded.vertexSize;
    size_t indexCount = encoded.indexCount;

    // Decode vertex buffer
    std::vector<uint8_t> interleavedVertices(vertexCount * stride);
    int vbResult = meshopt_decodeVertexBuffer(
        interleavedVertices.data(), vertexCount, stride,
        encoded.vertexData.data(), encoded.vertexData.size()
    );
    if (vbResult != 0) return {};

    // Decode index buffer
    std::vector<uint32_t> indices(indexCount);
    int ibResult = meshopt_decodeIndexBuffer(
        indices.data(), indexCount, sizeof(uint32_t),
        encoded.indexData.data(), encoded.indexData.size()
    );
    if (ibResult != 0) return {};

    // De-interleave into MeshData
    MeshData result;
    result.positions.resize(vertexCount * 3);
    if (hasNormals) result.normals.resize(vertexCount * 3);
    if (hasUVs) result.uvs.resize(vertexCount * 2);
    if (hasColors) result.colors.resize(vertexCount * 4);
    result.indices = std::move(indices);

    for (size_t v = 0; v < vertexCount; ++v) {
        const uint8_t* src = interleavedVertices.data() + v * stride;
        size_t off = 0;

        std::memcpy(&result.positions[v * 3], src + off, 12);
        off += 12;

        if (hasNormals) {
            std::memcpy(&result.normals[v * 3], src + off, 12);
            off += 12;
        }
        if (hasUVs) {
            std::memcpy(&result.uvs[v * 2], src + off, 8);
            off += 8;
        }
        if (hasColors) {
            std::memcpy(&result.colors[v * 4], src + off, 16);
            off += 16;
        }
    }

    return result;
#else
    (void)encoded; (void)hasNormals; (void)hasUVs; (void)hasColors;
    return {};
#endif
}

std::vector<uint8_t> encodeIndexBuffer(const std::vector<uint32_t>& indices, size_t vertexCount) {
#ifdef BROMESH_HAS_MESHOPTIMIZER
    if (indices.empty()) return {};

    size_t bound = meshopt_encodeIndexBufferBound(indices.size(), vertexCount);
    std::vector<uint8_t> encoded(bound);
    size_t size = meshopt_encodeIndexBuffer(
        encoded.data(), encoded.size(),
        indices.data(), indices.size()
    );
    encoded.resize(size);
    return encoded;
#else
    (void)indices; (void)vertexCount;
    return {};
#endif
}

std::vector<uint32_t> decodeIndexBuffer(const std::vector<uint8_t>& data, size_t indexCount) {
#ifdef BROMESH_HAS_MESHOPTIMIZER
    if (data.empty()) return {};

    std::vector<uint32_t> indices(indexCount);
    int result = meshopt_decodeIndexBuffer(
        indices.data(), indexCount, sizeof(uint32_t),
        data.data(), data.size()
    );
    if (result != 0) return {};
    return indices;
#else
    (void)data; (void)indexCount;
    return {};
#endif
}

} // namespace bromesh
