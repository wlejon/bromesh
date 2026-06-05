#pragma once

// Shared PLY header parsing + scalar reads, used by both the mesh loader
// (ply.cpp) and the Gaussian-splat loader (splat_ply.cpp). A 3DGS .ply is an
// ordinary PLY whose vertex element carries extra per-vertex float properties
// (f_dc_*, f_rest_*, opacity, scale_*, rot_*), so both loaders share the exact
// same header machinery and only differ in which properties they consume.
//
// Header-only with `inline` linkage: two TUs include it, so non-inline
// definitions would violate the ODR.

#include <cstdint>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

namespace bromesh::ply_detail {

enum class PlyFormat { ASCII, BinaryLE, BinaryBE };

struct PlyProp {
    std::string name;
    std::string type;      // float, uchar, uint, int, ...
    bool isList = false;
    std::string countType; // for list properties
    std::string valueType;
};

struct PlyHeader {
    PlyFormat format = PlyFormat::ASCII;
    size_t vertexCount = 0;
    size_t faceCount = 0;
    std::vector<PlyProp> vertexProps;
    std::vector<PlyProp> faceProps;
    size_t dataOffset = 0; // byte offset of the body (just past end_header\n)

    // Binary vertex layout (valid only when no vertex property is a list):
    // byte offset of each vertex property within one vertex record, and the
    // record stride. Empty for ASCII or when a vertex list property is present.
    std::vector<int> vertexOffsets;
    int vertexStride = 0;
};

inline int typeSize(const std::string& t) {
    if (t == "char" || t == "uchar" || t == "int8" || t == "uint8") return 1;
    if (t == "short" || t == "ushort" || t == "int16" || t == "uint16") return 2;
    if (t == "int" || t == "uint" || t == "int32" || t == "uint32" || t == "float" || t == "float32") return 4;
    if (t == "double" || t == "float64" || t == "int64" || t == "uint64") return 8;
    return 0;
}

inline float readFloatVal(const std::string& type, const uint8_t* data) {
    if (type == "float" || type == "float32") {
        float v; std::memcpy(&v, data, 4); return v;
    }
    if (type == "double" || type == "float64") {
        double v; std::memcpy(&v, data, 8); return static_cast<float>(v);
    }
    if (type == "uchar" || type == "uint8") return data[0] / 255.0f;
    if (type == "char" || type == "int8") return static_cast<int8_t>(data[0]) / 127.0f;
    if (type == "ushort" || type == "uint16") {
        uint16_t v; std::memcpy(&v, data, 2); return v / 65535.0f;
    }
    if (type == "int" || type == "int32") {
        int32_t v; std::memcpy(&v, data, 4); return static_cast<float>(v);
    }
    if (type == "uint" || type == "uint32") {
        uint32_t v; std::memcpy(&v, data, 4); return static_cast<float>(v);
    }
    return 0.0f;
}

inline uint32_t readUintVal(const std::string& type, const uint8_t* data) {
    if (type == "uchar" || type == "uint8") return data[0];
    if (type == "char" || type == "int8") return static_cast<uint8_t>(data[0]);
    if (type == "ushort" || type == "uint16") { uint16_t v; std::memcpy(&v, data, 2); return v; }
    if (type == "int" || type == "int32") { int32_t v; std::memcpy(&v, data, 4); return static_cast<uint32_t>(v); }
    if (type == "uint" || type == "uint32") { uint32_t v; std::memcpy(&v, data, 4); return v; }
    return 0;
}

/// Index of the vertex property with the given name, or -1.
inline int findProp(const std::vector<PlyProp>& props, const std::string& name) {
    for (int i = 0; i < static_cast<int>(props.size()); ++i)
        if (props[i].name == name) return i;
    return -1;
}

/// Parse a PLY header out of an in-memory file image. Fills format, element
/// counts, vertex/face property lists, and the body offset. For binary files
/// with no vertex list property, also fills vertexOffsets + vertexStride.
/// Returns false if no end_header marker is found.
inline bool parseHeader(const std::vector<uint8_t>& fileData, PlyHeader& out) {
    if (fileData.size() < 11) return false;

    size_t headerEnd = 0;
    std::string headerStr;
    for (size_t i = 0; i + 10 <= fileData.size(); ++i) {
        if (std::memcmp(&fileData[i], "end_header", 10) == 0) {
            size_t j = i + 10;
            while (j < fileData.size() && fileData[j] != '\n') j++;
            headerEnd = j + 1;
            headerStr.assign(reinterpret_cast<const char*>(fileData.data()), headerEnd);
            break;
        }
    }
    if (headerEnd == 0) return false;
    out.dataOffset = headerEnd;

    bool inVertex = false, inFace = false;
    std::istringstream hdr(headerStr);
    std::string line;
    while (std::getline(hdr, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        std::istringstream ls(line);
        std::string token;
        ls >> token;

        if (token == "format") {
            std::string fmt;
            ls >> fmt;
            if (fmt == "ascii") out.format = PlyFormat::ASCII;
            else if (fmt == "binary_little_endian") out.format = PlyFormat::BinaryLE;
            else if (fmt == "binary_big_endian") out.format = PlyFormat::BinaryBE;
        } else if (token == "element") {
            std::string ename;
            size_t ecount;
            ls >> ename >> ecount;
            if (ename == "vertex") { out.vertexCount = ecount; inVertex = true; inFace = false; }
            else if (ename == "face") { out.faceCount = ecount; inFace = true; inVertex = false; }
            else { inVertex = false; inFace = false; }
        } else if (token == "property") {
            std::string t1;
            ls >> t1;
            PlyProp prop;
            if (t1 == "list") {
                prop.isList = true;
                ls >> prop.countType >> prop.valueType >> prop.name;
            } else {
                prop.type = t1;
                ls >> prop.name;
            }
            if (inVertex) out.vertexProps.push_back(prop);
            else if (inFace) out.faceProps.push_back(prop);
        }
    }

    // Precompute binary vertex offsets/stride when the vertex record is fixed
    // size (no list properties — the 3DGS and ordinary mesh cases).
    bool vertexHasList = false;
    for (const auto& p : out.vertexProps) vertexHasList |= p.isList;
    if (out.format != PlyFormat::ASCII && !vertexHasList) {
        out.vertexOffsets.resize(out.vertexProps.size());
        int off = 0;
        for (size_t i = 0; i < out.vertexProps.size(); ++i) {
            out.vertexOffsets[i] = off;
            off += typeSize(out.vertexProps[i].type);
        }
        out.vertexStride = off;
    }
    return true;
}

} // namespace bromesh::ply_detail
