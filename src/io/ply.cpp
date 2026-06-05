#include "bromesh/io/ply.h"
#include "ply_common.h"

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

namespace bromesh {

bool savePLY(const MeshData& mesh, const std::string& path) {
    if (mesh.empty()) return false;

    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;

    size_t vertCount = mesh.vertexCount();
    size_t triCount = mesh.triangleCount();
    bool hasN = mesh.hasNormals();
    bool hasUV = mesh.hasUVs();
    bool hasC = mesh.hasColors();

    // Write ASCII header
    std::string header = "ply\nformat binary_little_endian 1.0\n";
    header += "element vertex " + std::to_string(vertCount) + "\n";
    header += "property float x\nproperty float y\nproperty float z\n";
    if (hasN)
        header += "property float nx\nproperty float ny\nproperty float nz\n";
    if (hasUV)
        header += "property float s\nproperty float t\n";
    if (hasC)
        header += "property uchar red\nproperty uchar green\nproperty uchar blue\nproperty uchar alpha\n";
    header += "element face " + std::to_string(triCount) + "\n";
    header += "property list uchar uint vertex_indices\n";
    header += "end_header\n";

    std::fwrite(header.c_str(), 1, header.size(), f);

    // Write vertex data
    for (size_t v = 0; v < vertCount; ++v) {
        float pos[3] = {
            mesh.positions[v * 3 + 0],
            mesh.positions[v * 3 + 1],
            mesh.positions[v * 3 + 2]
        };
        std::fwrite(pos, sizeof(float), 3, f);

        if (hasN) {
            float n[3] = {
                mesh.normals[v * 3 + 0],
                mesh.normals[v * 3 + 1],
                mesh.normals[v * 3 + 2]
            };
            std::fwrite(n, sizeof(float), 3, f);
        }
        if (hasUV) {
            float uv[2] = { mesh.uvs[v * 2 + 0], mesh.uvs[v * 2 + 1] };
            std::fwrite(uv, sizeof(float), 2, f);
        }
        if (hasC) {
            uint8_t col[4] = {
                static_cast<uint8_t>(std::min(255.0f, std::max(0.0f, mesh.colors[v * 4 + 0] * 255.0f))),
                static_cast<uint8_t>(std::min(255.0f, std::max(0.0f, mesh.colors[v * 4 + 1] * 255.0f))),
                static_cast<uint8_t>(std::min(255.0f, std::max(0.0f, mesh.colors[v * 4 + 2] * 255.0f))),
                static_cast<uint8_t>(std::min(255.0f, std::max(0.0f, mesh.colors[v * 4 + 3] * 255.0f)))
            };
            std::fwrite(col, 1, 4, f);
        }
    }

    // Write face data
    for (size_t t = 0; t < triCount; ++t) {
        uint8_t count = 3;
        std::fwrite(&count, 1, 1, f);
        uint32_t tri[3] = {
            mesh.indices[t * 3 + 0],
            mesh.indices[t * 3 + 1],
            mesh.indices[t * 3 + 2]
        };
        std::fwrite(tri, sizeof(uint32_t), 3, f);
    }

    std::fclose(f);
    return true;
}

// ---- PLY loading ----

using ply_detail::PlyFormat;
using ply_detail::PlyProp;
using ply_detail::typeSize;
using ply_detail::readFloatVal;
using ply_detail::readUintVal;

MeshData loadPLY(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return {};

    // Read entire file
    std::fseek(f, 0, SEEK_END);
    long fileSize = std::ftell(f);
    std::fseek(f, 0, SEEK_SET);
    std::vector<uint8_t> fileData(fileSize);
    if (static_cast<long>(std::fread(fileData.data(), 1, fileSize, f)) != fileSize) {
        std::fclose(f);
        return {};
    }
    std::fclose(f);

    // Parse header (shared with the Gaussian-splat loader)
    ply_detail::PlyHeader header;
    if (!ply_detail::parseHeader(fileData, header)) return {};

    const PlyFormat format = header.format;
    const size_t vertexCount = header.vertexCount;
    const size_t faceCount = header.faceCount;
    const std::vector<PlyProp>& vertexProps = header.vertexProps;
    const std::vector<PlyProp>& faceProps = header.faceProps;
    const size_t headerEnd = header.dataOffset;

    if (vertexCount == 0) return {};

    // Find property indices
    int xIdx = -1, yIdx = -1, zIdx = -1;
    int nxIdx = -1, nyIdx = -1, nzIdx = -1;
    int sIdx = -1, tIdx = -1;
    int rIdx = -1, gIdx = -1, bIdx = -1, aIdx = -1;

    for (int i = 0; i < static_cast<int>(vertexProps.size()); ++i) {
        const auto& p = vertexProps[i];
        if (p.name == "x") xIdx = i;
        else if (p.name == "y") yIdx = i;
        else if (p.name == "z") zIdx = i;
        else if (p.name == "nx") nxIdx = i;
        else if (p.name == "ny") nyIdx = i;
        else if (p.name == "nz") nzIdx = i;
        else if (p.name == "s" || p.name == "texture_u" || p.name == "u") sIdx = i;
        else if (p.name == "t" || p.name == "texture_v" || p.name == "v") tIdx = i;
        else if (p.name == "red") rIdx = i;
        else if (p.name == "green") gIdx = i;
        else if (p.name == "blue") bIdx = i;
        else if (p.name == "alpha") aIdx = i;
    }

    if (xIdx < 0 || yIdx < 0 || zIdx < 0) return {};

    bool hasN = nxIdx >= 0 && nyIdx >= 0 && nzIdx >= 0;
    bool hasUV = sIdx >= 0 && tIdx >= 0;
    bool hasC = rIdx >= 0 && gIdx >= 0 && bIdx >= 0;

    MeshData mesh;
    mesh.positions.resize(vertexCount * 3);
    if (hasN) mesh.normals.resize(vertexCount * 3);
    if (hasUV) mesh.uvs.resize(vertexCount * 2);
    if (hasC) mesh.colors.resize(vertexCount * 4, 1.0f);

    if (format == PlyFormat::BinaryLE || format == PlyFormat::BinaryBE) {
        // Binary loading (only LE supported fully; BE is rare)
        if (format == PlyFormat::BinaryBE) return {}; // Unsupported for now

        // Stride/offsets were precomputed by parseHeader (empty if a vertex
        // list property is present, which this loader does not support).
        if (header.vertexStride == 0) return {};
        const int vertexStride = header.vertexStride;
        const std::vector<int>& propOffsets = header.vertexOffsets;

        const uint8_t* data = fileData.data() + headerEnd;
        size_t remaining = fileData.size() - headerEnd;

        // Read vertices
        if (remaining < vertexCount * static_cast<size_t>(vertexStride)) return {};
        for (size_t v = 0; v < vertexCount; ++v) {
            const uint8_t* vdata = data + v * vertexStride;
            mesh.positions[v * 3 + 0] = readFloatVal(vertexProps[xIdx].type, vdata + propOffsets[xIdx]);
            mesh.positions[v * 3 + 1] = readFloatVal(vertexProps[yIdx].type, vdata + propOffsets[yIdx]);
            mesh.positions[v * 3 + 2] = readFloatVal(vertexProps[zIdx].type, vdata + propOffsets[zIdx]);

            if (hasN) {
                mesh.normals[v * 3 + 0] = readFloatVal(vertexProps[nxIdx].type, vdata + propOffsets[nxIdx]);
                mesh.normals[v * 3 + 1] = readFloatVal(vertexProps[nyIdx].type, vdata + propOffsets[nyIdx]);
                mesh.normals[v * 3 + 2] = readFloatVal(vertexProps[nzIdx].type, vdata + propOffsets[nzIdx]);
            }
            if (hasUV) {
                mesh.uvs[v * 2 + 0] = readFloatVal(vertexProps[sIdx].type, vdata + propOffsets[sIdx]);
                mesh.uvs[v * 2 + 1] = readFloatVal(vertexProps[tIdx].type, vdata + propOffsets[tIdx]);
            }
            if (hasC) {
                mesh.colors[v * 4 + 0] = readFloatVal(vertexProps[rIdx].type, vdata + propOffsets[rIdx]);
                mesh.colors[v * 4 + 1] = readFloatVal(vertexProps[gIdx].type, vdata + propOffsets[gIdx]);
                mesh.colors[v * 4 + 2] = readFloatVal(vertexProps[bIdx].type, vdata + propOffsets[bIdx]);
                if (aIdx >= 0)
                    mesh.colors[v * 4 + 3] = readFloatVal(vertexProps[aIdx].type, vdata + propOffsets[aIdx]);
            }
        }

        // Read faces
        const uint8_t* fdata = data + vertexCount * vertexStride;
        size_t faceDataRemaining = remaining - vertexCount * vertexStride;
        size_t foff = 0;

        // Find the face list property
        int faceListIdx = -1;
        for (int i = 0; i < static_cast<int>(faceProps.size()); ++i) {
            if (faceProps[i].isList) { faceListIdx = i; break; }
        }

        if (faceListIdx >= 0) {
            const auto& fp = faceProps[faceListIdx];
            int countSize = typeSize(fp.countType);
            int valSize = typeSize(fp.valueType);

            for (size_t fi = 0; fi < faceCount && foff < faceDataRemaining; ++fi) {
                uint32_t numVerts = readUintVal(fp.countType, fdata + foff);
                foff += countSize;

                if (numVerts >= 3) {
                    uint32_t first = readUintVal(fp.valueType, fdata + foff);
                    uint32_t prev = readUintVal(fp.valueType, fdata + foff + valSize);
                    for (uint32_t vi = 2; vi < numVerts; ++vi) {
                        uint32_t curr = readUintVal(fp.valueType, fdata + foff + vi * valSize);
                        mesh.indices.push_back(first);
                        mesh.indices.push_back(prev);
                        mesh.indices.push_back(curr);
                        prev = curr;
                    }
                }
                foff += numVerts * valSize;
            }
        }
    } else {
        // ASCII format
        std::string asciiData(reinterpret_cast<char*>(fileData.data() + headerEnd),
                              fileData.size() - headerEnd);
        std::istringstream ss(asciiData);

        for (size_t v = 0; v < vertexCount; ++v) {
            std::vector<float> vals(vertexProps.size());
            for (size_t i = 0; i < vertexProps.size(); ++i)
                ss >> vals[i];

            mesh.positions[v * 3 + 0] = vals[xIdx];
            mesh.positions[v * 3 + 1] = vals[yIdx];
            mesh.positions[v * 3 + 2] = vals[zIdx];

            if (hasN) {
                mesh.normals[v * 3 + 0] = vals[nxIdx];
                mesh.normals[v * 3 + 1] = vals[nyIdx];
                mesh.normals[v * 3 + 2] = vals[nzIdx];
            }
            if (hasUV) {
                mesh.uvs[v * 2 + 0] = vals[sIdx];
                mesh.uvs[v * 2 + 1] = vals[tIdx];
            }
            if (hasC) {
                mesh.colors[v * 4 + 0] = vals[rIdx] / 255.0f;
                mesh.colors[v * 4 + 1] = vals[gIdx] / 255.0f;
                mesh.colors[v * 4 + 2] = vals[bIdx] / 255.0f;
                if (aIdx >= 0)
                    mesh.colors[v * 4 + 3] = vals[aIdx] / 255.0f;
            }
        }

        for (size_t fi = 0; fi < faceCount; ++fi) {
            int numVerts;
            ss >> numVerts;
            std::vector<uint32_t> faceVerts(numVerts);
            for (int i = 0; i < numVerts; ++i)
                ss >> faceVerts[i];

            // Fan triangulation
            for (int i = 2; i < numVerts; ++i) {
                mesh.indices.push_back(faceVerts[0]);
                mesh.indices.push_back(faceVerts[i - 1]);
                mesh.indices.push_back(faceVerts[i]);
            }
        }
    }

    return mesh;
}

} // namespace bromesh
