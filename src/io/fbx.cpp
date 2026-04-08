#include "bromesh/io/fbx.h"

#ifdef BROMESH_HAS_OPENFBX
#include <ofbx.h>
#include <cstdio>
#include <cstring>
#include <memory>
#include <vector>
#include <algorithm>

namespace bromesh {

std::vector<MeshData> loadFBX(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return {};

    std::fseek(f, 0, SEEK_END);
    long fileSize = std::ftell(f);
    std::fseek(f, 0, SEEK_SET);

    std::vector<ofbx::u8> data(fileSize);
    if (static_cast<long>(std::fread(data.data(), 1, fileSize, f)) != fileSize) {
        std::fclose(f);
        return {};
    }
    std::fclose(f);

    ofbx::IScene* scene = ofbx::load(
        data.data(), static_cast<ofbx::usize>(fileSize),
        static_cast<ofbx::u16>(ofbx::LoadFlags::NONE));

    if (!scene) return {};

    // RAII guard for scene->destroy()
    struct SceneGuard {
        ofbx::IScene* s;
        ~SceneGuard() { if (s) s->destroy(); }
    } guard{scene};

    int meshCount = scene->getMeshCount();
    std::vector<MeshData> result;
    result.reserve(meshCount);

    for (int mi = 0; mi < meshCount; ++mi) {
        const ofbx::Mesh* fbxMesh = scene->getMesh(mi);
        if (!fbxMesh) continue;

        const ofbx::GeometryData& geoData = fbxMesh->getGeometryData();
        if (!geoData.hasVertices()) continue;

        auto positions = geoData.getPositions();
        auto normals = geoData.getNormals();
        auto uvs = geoData.getUVs(0);
        auto colors = geoData.getColors();

        bool hasNormals = normals.values != nullptr && normals.count > 0;
        bool hasUVs = uvs.values != nullptr && uvs.count > 0;
        bool hasColors = colors.values != nullptr && colors.count > 0;

        // Iterate partitions to triangulate polygons
        int partCount = geoData.getPartitionCount();
        if (partCount == 0) continue;

        MeshData mesh;

        for (int pi = 0; pi < partCount; ++pi) {
            auto partition = geoData.getPartition(pi);

            for (int poly = 0; poly < partition.polygon_count; ++poly) {
                const auto& p = partition.polygons[poly];
                int numVerts = p.vertex_count;

                // Triangulate polygon as a fan
                for (int vi = 2; vi < numVerts; ++vi) {
                    int indices[3] = {
                        p.from_vertex,
                        p.from_vertex + vi - 1,
                        p.from_vertex + vi
                    };

                    for (int k = 0; k < 3; ++k) {
                        int idx = indices[k];
                        uint32_t outIdx = static_cast<uint32_t>(mesh.positions.size() / 3);

                        // Position
                        int posIdx = positions.indices ? positions.indices[idx] : idx;
                        if (posIdx < 0) posIdx = -posIdx - 1; // FBX negative index convention
                        mesh.positions.push_back(positions.values[posIdx].x);
                        mesh.positions.push_back(positions.values[posIdx].y);
                        mesh.positions.push_back(positions.values[posIdx].z);

                        // Normal
                        if (hasNormals) {
                            int nIdx = normals.indices ? normals.indices[idx] : idx;
                            if (nIdx < 0) nIdx = -nIdx - 1;
                            if (nIdx < normals.count) {
                                mesh.normals.push_back(normals.values[nIdx].x);
                                mesh.normals.push_back(normals.values[nIdx].y);
                                mesh.normals.push_back(normals.values[nIdx].z);
                            } else {
                                mesh.normals.push_back(0); mesh.normals.push_back(1); mesh.normals.push_back(0);
                            }
                        }

                        // UV
                        if (hasUVs) {
                            int uIdx = uvs.indices ? uvs.indices[idx] : idx;
                            if (uIdx < 0) uIdx = -uIdx - 1;
                            if (uIdx < uvs.count) {
                                mesh.uvs.push_back(uvs.values[uIdx].x);
                                mesh.uvs.push_back(uvs.values[uIdx].y);
                            } else {
                                mesh.uvs.push_back(0); mesh.uvs.push_back(0);
                            }
                        }

                        // Colors
                        if (hasColors) {
                            int cIdx = colors.indices ? colors.indices[idx] : idx;
                            if (cIdx < 0) cIdx = -cIdx - 1;
                            if (cIdx < colors.count) {
                                mesh.colors.push_back(colors.values[cIdx].x);
                                mesh.colors.push_back(colors.values[cIdx].y);
                                mesh.colors.push_back(colors.values[cIdx].z);
                                mesh.colors.push_back(colors.values[cIdx].w);
                            } else {
                                mesh.colors.push_back(1); mesh.colors.push_back(1);
                                mesh.colors.push_back(1); mesh.colors.push_back(1);
                            }
                        }

                        mesh.indices.push_back(outIdx);
                    }
                }
            }
        }

        if (!mesh.empty()) {
            result.push_back(std::move(mesh));
        }
    }

    return result;
}

} // namespace bromesh

#else // !BROMESH_HAS_OPENFBX

namespace bromesh {

std::vector<MeshData> loadFBX(const std::string& /*path*/) {
    return {};
}

} // namespace bromesh

#endif
