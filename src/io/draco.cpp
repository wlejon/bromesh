#include "bromesh/io/draco.h"

#include "draco/compression/decode.h"
#include "draco/compression/encode.h"
#include "draco/mesh/triangle_soup_mesh_builder.h"

#include <cstring>
#include <memory>

namespace bromesh {

namespace {

const char* typeNameOf(draco::GeometryAttribute::Type t) {
    switch (t) {
        case draco::GeometryAttribute::POSITION: return "POSITION";
        case draco::GeometryAttribute::NORMAL: return "NORMAL";
        case draco::GeometryAttribute::COLOR: return "COLOR";
        case draco::GeometryAttribute::TEX_COORD: return "TEX_COORD";
        default: return "GENERIC";
    }
}

bool kindOf(draco::DataType dt, DracoAttribute::Kind& out) {
    switch (dt) {
        case draco::DT_FLOAT32: out = DracoAttribute::Kind::Float32; return true;
        case draco::DT_INT8: out = DracoAttribute::Kind::Int8; return true;
        case draco::DT_UINT8: out = DracoAttribute::Kind::Uint8; return true;
        case draco::DT_INT16: out = DracoAttribute::Kind::Int16; return true;
        case draco::DT_UINT16: out = DracoAttribute::Kind::Uint16; return true;
        case draco::DT_INT32: out = DracoAttribute::Kind::Int32; return true;
        case draco::DT_UINT32: out = DracoAttribute::Kind::Uint32; return true;
        default: return false;
    }
}

size_t kindSize(DracoAttribute::Kind k) {
    switch (k) {
        case DracoAttribute::Kind::Int8:
        case DracoAttribute::Kind::Uint8: return 1;
        case DracoAttribute::Kind::Int16:
        case DracoAttribute::Kind::Uint16: return 2;
        default: return 4;
    }
}

/// The attribute's values for every point, in point order, converted to
/// float. Returns false if the conversion is unsupported.
bool floatsOf(const draco::PointAttribute& att, uint32_t pointCount, int components,
              std::vector<float>& out) {
    out.resize(static_cast<size_t>(pointCount) * components);
    for (uint32_t p = 0; p < pointCount; ++p) {
        const draco::AttributeValueIndex avi = att.mapped_index(draco::PointIndex(p));
        if (!att.ConvertValue<float>(avi, static_cast<int8_t>(components),
                                     out.data() + static_cast<size_t>(p) * components)) {
            return false;
        }
    }
    return true;
}

}  // namespace

DracoDecoded decodeDraco(const uint8_t* data, size_t size) {
    DracoDecoded out;
    if (!data || size == 0) {
        out.error = "draco: empty buffer";
        return out;
    }

    draco::DecoderBuffer buffer;
    buffer.Init(reinterpret_cast<const char*>(data), size);

    auto typeStatus = draco::Decoder::GetEncodedGeometryType(&buffer);
    if (!typeStatus.ok()) {
        out.error = std::string("draco: not a draco buffer: ") +
                    typeStatus.status().error_msg();
        return out;
    }

    draco::Decoder decoder;
    std::unique_ptr<draco::PointCloud> cloud;
    const draco::Mesh* mesh = nullptr;
    if (typeStatus.value() == draco::TRIANGULAR_MESH) {
        auto status = decoder.DecodeMeshFromBuffer(&buffer);
        if (!status.ok()) {
            out.error = std::string("draco: mesh decode failed: ") +
                        status.status().error_msg();
            return out;
        }
        std::unique_ptr<draco::Mesh> m = std::move(status).value();
        mesh = m.get();
        cloud = std::move(m);
    } else if (typeStatus.value() == draco::POINT_CLOUD) {
        auto status = decoder.DecodePointCloudFromBuffer(&buffer);
        if (!status.ok()) {
            out.error = std::string("draco: point cloud decode failed: ") +
                        status.status().error_msg();
            return out;
        }
        cloud = std::move(status).value();
    } else {
        out.error = "draco: unknown geometry type";
        return out;
    }

    const uint32_t pointCount = cloud->num_points();

    // Every attribute, raw: the component type the decoder produced, one
    // value per POINT (draco's value array may be smaller and shared through
    // the point mapping, which a consumer should not have to know about).
    for (int i = 0; i < cloud->num_attributes(); ++i) {
        const draco::PointAttribute* att = cloud->attribute(i);
        DracoAttribute::Kind kind;
        if (!att || !kindOf(att->data_type(), kind)) continue;
        DracoAttribute& a = out.attributes.emplace_back();
        a.type = typeNameOf(att->attribute_type());
        a.uniqueId = att->unique_id();
        a.components = att->num_components();
        a.kind = kind;
        a.count = pointCount;
        const size_t valueBytes = kindSize(kind) * static_cast<size_t>(a.components);
        a.bytes.resize(valueBytes * pointCount);
        for (uint32_t p = 0; p < pointCount; ++p) {
            const draco::AttributeValueIndex avi =
                att->mapped_index(draco::PointIndex(p));
            std::memcpy(a.bytes.data() + valueBytes * p, att->GetAddress(avi),
                        valueBytes);
        }
    }

    // The standard streams, converted to MeshData's shape. First of each
    // semantic wins, which matches how the format is used — a second UV set
    // is still in `attributes`.
    const draco::PointAttribute* position =
        cloud->GetNamedAttribute(draco::GeometryAttribute::POSITION);
    if (!position) {
        out.error = "draco: no POSITION attribute";
        out.attributes.clear();
        return out;
    }
    if (!floatsOf(*position, pointCount, 3, out.mesh.positions)) {
        out.error = "draco: POSITION conversion failed";
        out.attributes.clear();
        return out;
    }

    if (const draco::PointAttribute* normal =
            cloud->GetNamedAttribute(draco::GeometryAttribute::NORMAL)) {
        if (normal->num_components() == 3)
            floatsOf(*normal, pointCount, 3, out.mesh.normals);
    }
    if (const draco::PointAttribute* uv =
            cloud->GetNamedAttribute(draco::GeometryAttribute::TEX_COORD)) {
        if (uv->num_components() >= 2) floatsOf(*uv, pointCount, 2, out.mesh.uvs);
    }
    if (const draco::PointAttribute* color =
            cloud->GetNamedAttribute(draco::GeometryAttribute::COLOR)) {
        const int comps = color->num_components();
        if (comps == 3 || comps == 4) {
            std::vector<float> raw;
            if (floatsOf(*color, pointCount, comps, raw)) {
                // An integer color converts to its numeric value; rgba wants
                // [0,1], so a quantized color normalizes by its type's max.
                float scale = 1.0f;
                switch (color->data_type()) {
                    case draco::DT_INT8: scale = 1.0f / 127.0f; break;
                    case draco::DT_UINT8: scale = 1.0f / 255.0f; break;
                    case draco::DT_INT16: scale = 1.0f / 32767.0f; break;
                    case draco::DT_UINT16: scale = 1.0f / 65535.0f; break;
                    default: break;
                }
                out.mesh.colors.resize(static_cast<size_t>(pointCount) * 4);
                for (uint32_t p = 0; p < pointCount; ++p) {
                    for (int c = 0; c < 4; ++c) {
                        out.mesh.colors[static_cast<size_t>(p) * 4 + c] =
                            c < comps ? raw[static_cast<size_t>(p) * comps + c] * scale
                                      : 1.0f;
                    }
                }
            }
        }
    }

    if (mesh) {
        out.mesh.indices.reserve(static_cast<size_t>(mesh->num_faces()) * 3);
        for (draco::FaceIndex f(0); f < mesh->num_faces(); ++f) {
            const draco::Mesh::Face& face = mesh->face(f);
            out.mesh.indices.push_back(face[0].value());
            out.mesh.indices.push_back(face[1].value());
            out.mesh.indices.push_back(face[2].value());
        }
    }

    return out;
}

std::vector<uint8_t> encodeDraco(const MeshData& mesh, const DracoEncodeOptions& options,
                                 std::string* error) {
    auto fail = [&](const char* why) {
        if (error) *error = why;
        return std::vector<uint8_t>{};
    };
    if (!mesh.validate()) return fail("draco: mesh fails validate()");
    if (mesh.positions.empty()) return fail("draco: no positions");
    if (mesh.indices.empty()) return fail("draco: no indices (only triangle meshes encode)");

    // Triangle soup in, connectivity rebuilt by the builder: each corner's
    // values are looked up through the index buffer, and the builder welds
    // identical corners back together.
    const uint32_t faceCount = static_cast<uint32_t>(mesh.indices.size() / 3);
    draco::TriangleSoupMeshBuilder builder;
    builder.Start(faceCount);

    const int posId = builder.AddAttribute(draco::GeometryAttribute::POSITION, 3,
                                           draco::DT_FLOAT32);
    const bool hasNormals = mesh.hasNormals();
    const bool hasUVs = mesh.hasUVs();
    const bool hasColors = mesh.hasColors();
    const int nrmId = hasNormals ? builder.AddAttribute(draco::GeometryAttribute::NORMAL,
                                                        3, draco::DT_FLOAT32)
                                 : -1;
    const int uvId = hasUVs ? builder.AddAttribute(draco::GeometryAttribute::TEX_COORD, 2,
                                                   draco::DT_FLOAT32)
                            : -1;
    const int colId = hasColors ? builder.AddAttribute(draco::GeometryAttribute::COLOR, 4,
                                                       draco::DT_FLOAT32)
                                : -1;

    for (uint32_t f = 0; f < faceCount; ++f) {
        const uint32_t i0 = mesh.indices[f * 3 + 0];
        const uint32_t i1 = mesh.indices[f * 3 + 1];
        const uint32_t i2 = mesh.indices[f * 3 + 2];
        const draco::FaceIndex face(f);
        builder.SetAttributeValuesForFace(posId, face, &mesh.positions[i0 * 3],
                                          &mesh.positions[i1 * 3], &mesh.positions[i2 * 3]);
        if (hasNormals) {
            builder.SetAttributeValuesForFace(nrmId, face, &mesh.normals[i0 * 3],
                                              &mesh.normals[i1 * 3], &mesh.normals[i2 * 3]);
        }
        if (hasUVs) {
            builder.SetAttributeValuesForFace(uvId, face, &mesh.uvs[i0 * 2],
                                              &mesh.uvs[i1 * 2], &mesh.uvs[i2 * 2]);
        }
        if (hasColors) {
            builder.SetAttributeValuesForFace(colId, face, &mesh.colors[i0 * 4],
                                              &mesh.colors[i1 * 4], &mesh.colors[i2 * 4]);
        }
    }

    std::unique_ptr<draco::Mesh> dracoMesh = builder.Finalize();
    if (!dracoMesh) return fail("draco: mesh build failed");

    draco::Encoder encoder;
    encoder.SetAttributeQuantization(draco::GeometryAttribute::POSITION,
                                     options.positionBits);
    encoder.SetAttributeQuantization(draco::GeometryAttribute::NORMAL, options.normalBits);
    encoder.SetAttributeQuantization(draco::GeometryAttribute::TEX_COORD, options.uvBits);
    encoder.SetAttributeQuantization(draco::GeometryAttribute::COLOR, options.colorBits);
    encoder.SetSpeedOptions(options.speed, options.speed);

    draco::EncoderBuffer buffer;
    const draco::Status status = encoder.EncodeMeshToBuffer(*dracoMesh, &buffer);
    if (!status.ok()) {
        if (error) *error = std::string("draco: encode failed: ") + status.error_msg();
        return {};
    }
    const uint8_t* bytes = reinterpret_cast<const uint8_t*>(buffer.data());
    return std::vector<uint8_t>(bytes, bytes + buffer.size());
}

}  // namespace bromesh
