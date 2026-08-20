#include "test_framework.h"
#include <algorithm>
#include <cmath>
#include <cstring>

#if BROMESH_HAS_DRACO

#include "bromesh/io/draco.h"
#include "bromesh/primitives/primitives.h"

using bromesh::MeshData;

// A quantized round trip cannot promise exact floats; positions at 14 bits
// over a unit-scale mesh land within ~1e-3.
static bool near(float a, float b, float tol) { return std::fabs(a - b) <= tol; }

TEST(draco_round_trip_positions_and_indices) {
    MeshData box = bromesh::box(1.0f, 0.5f, 0.25f);
    ASSERT(box.validate(), "box validates");

    std::string error;
    std::vector<uint8_t> bytes = bromesh::encodeDraco(box, bromesh::DracoEncodeOptions(), &error);
    ASSERT(!bytes.empty(), ("encode produced bytes: " + error).c_str());

    bromesh::DracoDecoded decoded = bromesh::decodeDraco(bytes.data(), bytes.size());
    ASSERT(decoded.ok(), decoded.error.c_str());
    ASSERT(decoded.mesh.validate(), "decoded mesh validates");
    ASSERT(decoded.mesh.triangleCount() == box.triangleCount(),
           "triangle count survives the trip");

    // The soup builder rewelds, so vertex ORDER changes; the geometry must
    // not. Compare the axis-aligned extents, which quantization only nudges.
    float lo[3] = {1e9f, 1e9f, 1e9f}, hi[3] = {-1e9f, -1e9f, -1e9f};
    float dlo[3] = {1e9f, 1e9f, 1e9f}, dhi[3] = {-1e9f, -1e9f, -1e9f};
    for (size_t i = 0; i < box.positions.size(); i += 3)
        for (int c = 0; c < 3; ++c) {
            lo[c] = std::min(lo[c], box.positions[i + c]);
            hi[c] = std::max(hi[c], box.positions[i + c]);
        }
    for (size_t i = 0; i < decoded.mesh.positions.size(); i += 3)
        for (int c = 0; c < 3; ++c) {
            dlo[c] = std::min(dlo[c], decoded.mesh.positions[i + c]);
            dhi[c] = std::max(dhi[c], decoded.mesh.positions[i + c]);
        }
    for (int c = 0; c < 3; ++c) {
        ASSERT(near(lo[c], dlo[c], 1e-3f), "min extent survives");
        ASSERT(near(hi[c], dhi[c], 1e-3f), "max extent survives");
    }
}

TEST(draco_attributes_surface_raw) {
    MeshData sphere = bromesh::sphere(1.0f, 12, 8);
    ASSERT(sphere.hasNormals(), "sphere has normals");
    ASSERT(sphere.hasUVs(), "sphere has uvs");

    std::vector<uint8_t> bytes = bromesh::encodeDraco(sphere);
    ASSERT(!bytes.empty(), "encode produced bytes");
    bromesh::DracoDecoded decoded = bromesh::decodeDraco(bytes.data(), bytes.size());
    ASSERT(decoded.ok(), decoded.error.c_str());

    ASSERT(decoded.mesh.hasNormals(), "normals decode into the standard stream");
    ASSERT(decoded.mesh.hasUVs(), "uvs decode into the standard stream");

    bool sawPosition = false, sawNormal = false, sawUv = false;
    for (const bromesh::DracoAttribute& a : decoded.attributes) {
        if (a.type == "POSITION") sawPosition = true;
        if (a.type == "NORMAL") sawNormal = true;
        if (a.type == "TEX_COORD") sawUv = true;
        ASSERT(a.count == decoded.mesh.vertexCount(), "attribute count is vertex count");
        ASSERT(!a.bytes.empty(), "attribute carries bytes");
    }
    ASSERT(sawPosition && sawNormal && sawUv, "all three semantics in the raw list");

    // Unit normals must survive quantization as unit normals.
    for (size_t i = 0; i + 2 < decoded.mesh.normals.size(); i += 3) {
        const float len = std::sqrt(decoded.mesh.normals[i] * decoded.mesh.normals[i] +
                                    decoded.mesh.normals[i + 1] * decoded.mesh.normals[i + 1] +
                                    decoded.mesh.normals[i + 2] * decoded.mesh.normals[i + 2]);
        ASSERT(near(len, 1.0f, 0.05f), "normal stays unit length");
    }
}

TEST(draco_decode_rejects_garbage) {
    const uint8_t junk[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    bromesh::DracoDecoded decoded = bromesh::decodeDraco(junk, sizeof junk);
    ASSERT(!decoded.ok(), "garbage does not decode");
    ASSERT(bromesh::decodeDraco(nullptr, 0).ok() == false, "empty does not decode");
}

#endif  // BROMESH_HAS_DRACO
