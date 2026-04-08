#include "bromesh/io/gltf.h"

#if BROMESH_HAS_GLTF
#include "tiny_gltf.h"
#endif

namespace bromesh {

GltfScene loadGLTF(const std::string& path) {
#if BROMESH_HAS_GLTF
    tinygltf::Model model;
    tinygltf::TinyGLTF loader;
    std::string err, warn;

    bool ok = false;
    if (path.size() >= 4 && path.substr(path.size() - 4) == ".glb") {
        ok = loader.LoadBinaryFromFile(&model, &err, &warn, path);
    } else {
        ok = loader.LoadASCIIFromFile(&model, &err, &warn, path);
    }

    if (!ok) return {};

    GltfScene scene;

    for (const auto& mesh : model.meshes) {
        for (const auto& prim : mesh.primitives) {
            if (prim.mode != TINYGLTF_MODE_TRIANGLES) continue;

            MeshData md;

            // Indices
            if (prim.indices >= 0) {
                const auto& accessor = model.accessors[prim.indices];
                const auto& bufView = model.bufferViews[accessor.bufferView];
                const auto& buf = model.buffers[bufView.buffer];
                const uint8_t* base = buf.data.data() + bufView.byteOffset + accessor.byteOffset;

                md.indices.resize(accessor.count);
                if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT) {
                    const auto* src = reinterpret_cast<const uint16_t*>(base);
                    for (size_t i = 0; i < accessor.count; ++i)
                        md.indices[i] = src[i];
                } else if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT) {
                    const auto* src = reinterpret_cast<const uint32_t*>(base);
                    for (size_t i = 0; i < accessor.count; ++i)
                        md.indices[i] = src[i];
                }
            }

            // Positions
            {
                auto it = prim.attributes.find("POSITION");
                if (it != prim.attributes.end()) {
                    const auto& accessor = model.accessors[it->second];
                    const auto& bufView = model.bufferViews[accessor.bufferView];
                    const auto& buf = model.buffers[bufView.buffer];
                    const auto* src = reinterpret_cast<const float*>(
                        buf.data.data() + bufView.byteOffset + accessor.byteOffset);
                    md.positions.assign(src, src + accessor.count * 3);
                }
            }

            // Normals
            {
                auto it = prim.attributes.find("NORMAL");
                if (it != prim.attributes.end()) {
                    const auto& accessor = model.accessors[it->second];
                    const auto& bufView = model.bufferViews[accessor.bufferView];
                    const auto& buf = model.buffers[bufView.buffer];
                    const auto* src = reinterpret_cast<const float*>(
                        buf.data.data() + bufView.byteOffset + accessor.byteOffset);
                    md.normals.assign(src, src + accessor.count * 3);
                }
            }

            // UVs
            {
                auto it = prim.attributes.find("TEXCOORD_0");
                if (it != prim.attributes.end()) {
                    const auto& accessor = model.accessors[it->second];
                    const auto& bufView = model.bufferViews[accessor.bufferView];
                    const auto& buf = model.buffers[bufView.buffer];
                    const auto* src = reinterpret_cast<const float*>(
                        buf.data.data() + bufView.byteOffset + accessor.byteOffset);
                    md.uvs.assign(src, src + accessor.count * 2);
                }
            }

            scene.meshes.push_back(std::move(md));
        }
    }

    return scene;
#else
    (void)path;
    return {};
#endif
}

bool saveGLTF(const MeshData& mesh, const std::string& path) {
#if BROMESH_HAS_GLTF
    if (mesh.empty()) return false;

    tinygltf::Model model;
    tinygltf::Scene gltfScene;
    tinygltf::Mesh gltfMesh;
    tinygltf::Primitive prim;
    tinygltf::Node node;

    model.asset.version = "2.0";
    model.asset.generator = "bromesh";

    // Build a single buffer with positions, normals, and indices
    tinygltf::Buffer buffer;

    auto appendData = [&](const void* data, size_t bytes) -> size_t {
        size_t offset = buffer.data.size();
        buffer.data.resize(offset + bytes);
        memcpy(buffer.data.data() + offset, data, bytes);
        return offset;
    };

    // Positions
    size_t posOffset = appendData(mesh.positions.data(),
                                  mesh.positions.size() * sizeof(float));
    // Normals
    size_t normOffset = 0;
    if (mesh.hasNormals()) {
        normOffset = appendData(mesh.normals.data(),
                                mesh.normals.size() * sizeof(float));
    }
    // UVs
    size_t uvOffset = 0;
    if (mesh.hasUVs()) {
        uvOffset = appendData(mesh.uvs.data(),
                              mesh.uvs.size() * sizeof(float));
    }
    // Indices
    size_t idxOffset = appendData(mesh.indices.data(),
                                  mesh.indices.size() * sizeof(uint32_t));

    model.buffers.push_back(std::move(buffer));

    int bufViewIdx = 0;

    // Position buffer view
    {
        tinygltf::BufferView bv;
        bv.buffer = 0;
        bv.byteOffset = posOffset;
        bv.byteLength = mesh.positions.size() * sizeof(float);
        bv.target = TINYGLTF_TARGET_ARRAY_BUFFER;
        model.bufferViews.push_back(bv);
    }
    int posBvIdx = bufViewIdx++;

    // Position accessor
    {
        tinygltf::Accessor acc;
        acc.bufferView = posBvIdx;
        acc.byteOffset = 0;
        acc.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
        acc.type = TINYGLTF_TYPE_VEC3;
        acc.count = mesh.vertexCount();
        model.accessors.push_back(acc);
    }
    prim.attributes["POSITION"] = 0;

    int accIdx = 1;

    // Normals
    if (mesh.hasNormals()) {
        tinygltf::BufferView bv;
        bv.buffer = 0;
        bv.byteOffset = normOffset;
        bv.byteLength = mesh.normals.size() * sizeof(float);
        bv.target = TINYGLTF_TARGET_ARRAY_BUFFER;
        model.bufferViews.push_back(bv);
        int normBvIdx = bufViewIdx++;

        tinygltf::Accessor acc;
        acc.bufferView = normBvIdx;
        acc.byteOffset = 0;
        acc.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
        acc.type = TINYGLTF_TYPE_VEC3;
        acc.count = mesh.vertexCount();
        model.accessors.push_back(acc);
        prim.attributes["NORMAL"] = accIdx++;
    }

    // UVs
    if (mesh.hasUVs()) {
        tinygltf::BufferView bv;
        bv.buffer = 0;
        bv.byteOffset = uvOffset;
        bv.byteLength = mesh.uvs.size() * sizeof(float);
        bv.target = TINYGLTF_TARGET_ARRAY_BUFFER;
        model.bufferViews.push_back(bv);
        int uvBvIdx = bufViewIdx++;

        tinygltf::Accessor acc;
        acc.bufferView = uvBvIdx;
        acc.byteOffset = 0;
        acc.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
        acc.type = TINYGLTF_TYPE_VEC2;
        acc.count = mesh.vertexCount();
        model.accessors.push_back(acc);
        prim.attributes["TEXCOORD_0"] = accIdx++;
    }

    // Indices
    {
        tinygltf::BufferView bv;
        bv.buffer = 0;
        bv.byteOffset = idxOffset;
        bv.byteLength = mesh.indices.size() * sizeof(uint32_t);
        bv.target = TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER;
        model.bufferViews.push_back(bv);
        int idxBvIdx = bufViewIdx++;

        tinygltf::Accessor acc;
        acc.bufferView = idxBvIdx;
        acc.byteOffset = 0;
        acc.componentType = TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT;
        acc.type = TINYGLTF_TYPE_SCALAR;
        acc.count = mesh.indices.size();
        model.accessors.push_back(acc);
        prim.indices = accIdx++;
    }

    prim.mode = TINYGLTF_MODE_TRIANGLES;
    gltfMesh.primitives.push_back(prim);
    model.meshes.push_back(gltfMesh);

    node.mesh = 0;
    model.nodes.push_back(node);
    gltfScene.nodes.push_back(0);
    model.scenes.push_back(gltfScene);
    model.defaultScene = 0;

    tinygltf::TinyGLTF writer;
    bool binary = path.size() >= 4 && path.substr(path.size() - 4) == ".glb";
    return writer.WriteGltfSceneToFile(&model, path,
                                       true,  // embedImages
                                       true,  // embedBuffers
                                       true,  // prettyPrint
                                       binary);
#else
    (void)mesh; (void)path;
    return false;
#endif
}

} // namespace bromesh
