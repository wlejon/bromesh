#pragma once

#include "bromesh/mesh_data.h"

#include <cstdint>
#include <string>
#include <vector>

namespace bromesh {

/// One decoded Draco attribute, raw. The standard streams land in
/// DracoDecoded::mesh already converted; this list carries EVERY attribute
/// with its original component type intact, which is what a glTF consumer
/// needs — KHR_draco_mesh_compression maps glTF accessors to Draco attributes
/// by unique id, and skinning data (JOINTS_0 as uint8/uint16) must not be
/// forced through float.
struct DracoAttribute {
    /// Draco's semantic: "POSITION", "NORMAL", "COLOR", "TEX_COORD" or
    /// "GENERIC".
    std::string type;
    /// The id KHR_draco_mesh_compression keys on.
    uint32_t uniqueId = 0;
    /// Components per vertex (e.g. 3 for a position).
    int components = 0;
    /// The component type the decoder produced.
    enum class Kind : uint8_t { Float32, Int8, Uint8, Int16, Uint16, Int32, Uint32 };
    Kind kind = Kind::Float32;
    /// count * components values of `kind`, tightly packed.
    std::vector<uint8_t> bytes;
    /// Vertex count (same for every attribute of one mesh).
    uint32_t count = 0;
};

/// Result of decoding a Draco buffer. A point cloud decodes with empty
/// indices. `error` is non-empty exactly when the decode failed, in which
/// case everything else is empty.
struct DracoDecoded {
    MeshData mesh;                          ///< standard streams, converted
    std::vector<DracoAttribute> attributes; ///< every attribute, raw
    std::string error;
    bool ok() const { return error.empty(); }
};

/// Decode a Draco-compressed mesh or point cloud from memory — a .drc file's
/// bytes, or the buffer view KHR_draco_mesh_compression points at. The
/// standard attributes fill DracoDecoded::mesh (colors expand to rgba;
/// a 16-bit UV dequantizes to float), and every attribute — standard and
/// generic alike — also appears raw in `attributes`.
DracoDecoded decodeDraco(const uint8_t* data, size_t size);

/// Quantization/speed knobs for encodeDraco. Bits follow Draco's usual
/// ranges (1..30); speed trades compression for time, 0 = smallest,
/// 10 = fastest.
struct DracoEncodeOptions {
    int positionBits = 14;
    int normalBits = 10;
    int uvBits = 12;
    int colorBits = 8;
    int speed = 7;
};

/// Encode an indexed mesh to a Draco buffer (.drc bytes). Positions are
/// required; normals/uvs/colors encode when present. Returns empty and sets
/// `error` (when given) on failure.
std::vector<uint8_t> encodeDraco(const MeshData& mesh,
                                 const DracoEncodeOptions& options = DracoEncodeOptions(),
                                 std::string* error = nullptr);

} // namespace bromesh
