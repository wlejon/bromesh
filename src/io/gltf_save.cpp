#include "gltf_common.h"

#if BROMESH_HAS_GLTF
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#pragma clang diagnostic ignored "-Wmissing-field-initializers"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif

#define STB_IMAGE_WRITE_STATIC
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
#endif

#include <filesystem>

namespace bromesh {

#if !BROMESH_HAS_GLTF
bool saveGLTF(const GltfScene&, const std::string&) { return false; }
bool saveGLTF(const MeshData&, const std::string&) { return false; }
bool saveGLTF(const MeshData&, const SkinData*, const Skeleton*,
              const std::vector<Animation>&, const std::string&) { return false; }
bool saveGLTF(const MeshData&,
              std::optional<std::reference_wrapper<const SkinData>>,
              std::optional<std::reference_wrapper<const Skeleton>>,
              const std::vector<Animation>&, const std::string&) { return false; }
#else

struct SaveCtx {
    tinygltf::Model model;
    tinygltf::Buffer buffer;

    int appendView(const void* data, size_t bytes, int target = 0) {
        while (buffer.data.size() % 4 != 0) {
            buffer.data.push_back(0);
        }
        size_t offset = buffer.data.size();
        buffer.data.resize(offset + bytes);
        std::memcpy(buffer.data.data() + offset, data, bytes);
        tinygltf::BufferView bv;
        bv.buffer = 0;
        bv.byteOffset = offset;
        bv.byteLength = bytes;
        if (target) bv.target = target;
        model.bufferViews.push_back(bv);
        return static_cast<int>(model.bufferViews.size()) - 1;
    }

    int addAccessor(int view, int componentType, int type, size_t count,
                    const std::vector<double>& minVals = {},
                    const std::vector<double>& maxVals = {}) {
        tinygltf::Accessor a;
        a.bufferView = view;
        a.byteOffset = 0;
        a.componentType = componentType;
        a.type = type;
        a.count = count;
        a.minValues = minVals;
        a.maxValues = maxVals;
        model.accessors.push_back(a);
        return static_cast<int>(model.accessors.size()) - 1;
    }
};

bool saveGLTF(const GltfScene& scene, const std::string& path) {
    if (scene.meshes.empty() && scene.skeletons.empty()) return false;
    for (const auto& m : scene.meshes) {
        if (m.empty()) return false;
    }

    std::error_code ec;
    std::filesystem::path fsPath(path);
    if (fsPath.has_parent_path()) {
        std::filesystem::create_directories(fsPath.parent_path(), ec);
    }

    SaveCtx ctx;
    ctx.model.asset.version = "2.0";
    ctx.model.asset.generator = "bromesh";

    // --- Images and textures ---
    for (const auto& img : scene.images) {
        tinygltf::Image gImg;
        gImg.name = img.name;
        std::string mime = img.mimeType.empty() ? "image/png" : img.mimeType;

        bool isPng = (img.data.size() >= 8 &&
                      img.data[0] == 0x89 && img.data[1] == 'P' &&
                      img.data[2] == 'N' && img.data[3] == 'G');
        bool isJpg = (img.data.size() >= 3 &&
                      img.data[0] == 0xFF && img.data[1] == 0xD8 && img.data[2] == 0xFF);

        if (isPng || isJpg) {
            int v = ctx.appendView(img.data.data(), img.data.size(), 0);
            gImg.bufferView = v;
            gImg.mimeType = isPng ? "image/png" : "image/jpeg";
        } else if (img.width > 0 && img.height > 0 && !img.data.empty()) {
            int len = 0;
            unsigned char* pngMem = stbi_write_png_to_mem(
                img.data.data(), img.width * 4, img.width, img.height, 4, &len);
            if (pngMem && len > 0) {
                int v = ctx.appendView(pngMem, static_cast<size_t>(len), 0);
                free(pngMem);
                gImg.bufferView = v;
                gImg.mimeType = "image/png";
            } else {
                int v = ctx.appendView(img.data.data(), img.data.size(), 0);
                gImg.bufferView = v;
                gImg.mimeType = mime;
            }
        } else if (!img.data.empty()) {
            int v = ctx.appendView(img.data.data(), img.data.size(), 0);
            gImg.bufferView = v;
            gImg.mimeType = mime;
        }
        ctx.model.images.push_back(gImg);

        tinygltf::Texture gTex;
        gTex.source = static_cast<int>(ctx.model.images.size()) - 1;
        ctx.model.textures.push_back(gTex);
    }

    // --- Materials ---
    for (const auto& mat : scene.materials) {
        tinygltf::Material gMat;
        gMat.name = mat.name;
        gMat.pbrMetallicRoughness.baseColorFactor = {
            static_cast<double>(mat.baseColorFactor[0]),
            static_cast<double>(mat.baseColorFactor[1]),
            static_cast<double>(mat.baseColorFactor[2]),
            static_cast<double>(mat.baseColorFactor[3])
        };
        gMat.pbrMetallicRoughness.metallicFactor = static_cast<double>(mat.metallicFactor);
        gMat.pbrMetallicRoughness.roughnessFactor = static_cast<double>(mat.roughnessFactor);
        gMat.emissiveFactor = {
            static_cast<double>(mat.emissiveFactor[0]),
            static_cast<double>(mat.emissiveFactor[1]),
            static_cast<double>(mat.emissiveFactor[2])
        };

        if (mat.baseColorTexture >= 0 && mat.baseColorTexture < static_cast<int>(ctx.model.textures.size())) {
            gMat.pbrMetallicRoughness.baseColorTexture.index = mat.baseColorTexture;
        }
        if (mat.metallicRoughnessTexture >= 0 && mat.metallicRoughnessTexture < static_cast<int>(ctx.model.textures.size())) {
            gMat.pbrMetallicRoughness.metallicRoughnessTexture.index = mat.metallicRoughnessTexture;
        }
        if (mat.normalTexture >= 0 && mat.normalTexture < static_cast<int>(ctx.model.textures.size())) {
            gMat.normalTexture.index = mat.normalTexture;
        }
        if (mat.occlusionTexture >= 0 && mat.occlusionTexture < static_cast<int>(ctx.model.textures.size())) {
            gMat.occlusionTexture.index = mat.occlusionTexture;
        }
        if (mat.emissiveTexture >= 0 && mat.emissiveTexture < static_cast<int>(ctx.model.textures.size())) {
            gMat.emissiveTexture.index = mat.emissiveTexture;
        }
        ctx.model.materials.push_back(gMat);
    }

    // --- Skeletons and skins ---
    std::vector<int> skeletonFirstBoneNode;
    std::vector<int> skeletonSkinIndex;
    skeletonFirstBoneNode.reserve(scene.skeletons.size());
    skeletonSkinIndex.reserve(scene.skeletons.size());

    for (const auto& sk : scene.skeletons) {
        int firstBoneNode = static_cast<int>(ctx.model.nodes.size());
        skeletonFirstBoneNode.push_back(firstBoneNode);

        for (const auto& b : sk.bones) {
            tinygltf::Node n;
            n.name = b.name;
            n.translation = { b.localT[0], b.localT[1], b.localT[2] };
            n.rotation    = { b.localR[0], b.localR[1], b.localR[2], b.localR[3] };
            n.scale       = { b.localS[0], b.localS[1], b.localS[2] };
            ctx.model.nodes.push_back(n);
        }

        for (size_t b = 0; b < sk.bones.size(); ++b) {
            int p = sk.bones[b].parent;
            if (p >= 0 && p < static_cast<int>(sk.bones.size())) {
                ctx.model.nodes[firstBoneNode + p].children.push_back(firstBoneNode + static_cast<int>(b));
            }
        }

        std::vector<float> ibm(sk.bones.size() * 16);
        for (size_t b = 0; b < sk.bones.size(); ++b) {
            std::memcpy(&ibm[b * 16], sk.bones[b].inverseBind, 16 * sizeof(float));
        }
        int vIbm = ctx.appendView(ibm.data(), ibm.size() * sizeof(float));
        int aIbm = ctx.addAccessor(vIbm, TINYGLTF_COMPONENT_TYPE_FLOAT,
                                    TINYGLTF_TYPE_MAT4, sk.bones.size());

        tinygltf::Skin gskin;
        gskin.inverseBindMatrices = aIbm;
        for (size_t b = 0; b < sk.bones.size(); ++b) {
            gskin.joints.push_back(firstBoneNode + static_cast<int>(b));
        }
        ctx.model.skins.push_back(gskin);
        skeletonSkinIndex.push_back(static_cast<int>(ctx.model.skins.size()) - 1);
    }

    // --- Meshes ---
    std::vector<int> meshNodeIndices;
    meshNodeIndices.reserve(scene.meshes.size());

    for (size_t mi = 0; mi < scene.meshes.size(); ++mi) {
        const auto& mesh = scene.meshes[mi];
        tinygltf::Mesh gltfMesh;
        tinygltf::Primitive prim;
        prim.mode = TINYGLTF_MODE_TRIANGLES;

        // Positions
        {
            float mn[3] = { mesh.positions[0], mesh.positions[1], mesh.positions[2] };
            float mx[3] = { mn[0], mn[1], mn[2] };
            for (size_t i = 0; i < mesh.positions.size(); i += 3) {
                for (int k = 0; k < 3; ++k) {
                    mn[k] = std::min(mn[k], mesh.positions[i + k]);
                    mx[k] = std::max(mx[k], mesh.positions[i + k]);
                }
            }
            int v = ctx.appendView(mesh.positions.data(),
                                   mesh.positions.size() * sizeof(float),
                                   TINYGLTF_TARGET_ARRAY_BUFFER);
            int a = ctx.addAccessor(v, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_VEC3,
                                    mesh.vertexCount(),
                                    {(double)mn[0], (double)mn[1], (double)mn[2]},
                                    {(double)mx[0], (double)mx[1], (double)mx[2]});
            prim.attributes["POSITION"] = a;
        }

        // Normals
        if (mesh.hasNormals()) {
            int v = ctx.appendView(mesh.normals.data(),
                                   mesh.normals.size() * sizeof(float),
                                   TINYGLTF_TARGET_ARRAY_BUFFER);
            int a = ctx.addAccessor(v, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_VEC3,
                                    mesh.vertexCount());
            prim.attributes["NORMAL"] = a;
        }

        // UVs
        if (mesh.hasUVs()) {
            int v = ctx.appendView(mesh.uvs.data(),
                                   mesh.uvs.size() * sizeof(float),
                                   TINYGLTF_TARGET_ARRAY_BUFFER);
            int a = ctx.addAccessor(v, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_VEC2,
                                    mesh.vertexCount());
            prim.attributes["TEXCOORD_0"] = a;
        }

        // Colors
        if (mesh.hasColors()) {
            int v = ctx.appendView(mesh.colors.data(),
                                   mesh.colors.size() * sizeof(float),
                                   TINYGLTF_TARGET_ARRAY_BUFFER);
            int a = ctx.addAccessor(v, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_VEC4,
                                    mesh.vertexCount());
            prim.attributes["COLOR_0"] = a;
        }

        // Tangents
        if (mesh.hasTangents()) {
            int v = ctx.appendView(mesh.tangents.data(),
                                   mesh.tangents.size() * sizeof(float),
                                   TINYGLTF_TARGET_ARRAY_BUFFER);
            int a = ctx.addAccessor(v, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_VEC4,
                                    mesh.vertexCount());
            prim.attributes["TANGENT"] = a;
        }

        // Skinning
        int skIdx = -1;
        if (mi < scene.meshSkeleton.size()) {
            skIdx = scene.meshSkeleton[mi];
        } else if (scene.skeletons.size() == 1) {
            skIdx = 0;
        }

        bool skinned = (skIdx >= 0 && skIdx < static_cast<int>(scene.skeletons.size()) &&
                        mi < scene.skins.size() &&
                        scene.skins[mi].boneIndices.size() == mesh.vertexCount() * 4 &&
                        scene.skins[mi].boneWeights.size() == mesh.vertexCount() * 4);

        if (skinned) {
            std::vector<uint16_t> joints(scene.skins[mi].boneIndices.size());
            for (size_t j = 0; j < joints.size(); ++j)
                joints[j] = static_cast<uint16_t>(scene.skins[mi].boneIndices[j]);
            int v = ctx.appendView(joints.data(),
                                   joints.size() * sizeof(uint16_t),
                                   TINYGLTF_TARGET_ARRAY_BUFFER);
            int a = ctx.addAccessor(v, TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT,
                                    TINYGLTF_TYPE_VEC4, mesh.vertexCount());
            prim.attributes["JOINTS_0"] = a;

            int v2 = ctx.appendView(scene.skins[mi].boneWeights.data(),
                                    scene.skins[mi].boneWeights.size() * sizeof(float),
                                    TINYGLTF_TARGET_ARRAY_BUFFER);
            int a2 = ctx.addAccessor(v2, TINYGLTF_COMPONENT_TYPE_FLOAT,
                                     TINYGLTF_TYPE_VEC4, mesh.vertexCount());
            prim.attributes["WEIGHTS_0"] = a2;
        }

        // Material assignment
        if (mi < scene.meshMaterial.size() && scene.meshMaterial[mi] >= 0 &&
            scene.meshMaterial[mi] < static_cast<int>(ctx.model.materials.size())) {
            prim.material = scene.meshMaterial[mi];
        } else if (scene.meshMaterial.empty() && scene.materials.size() == 1 && scene.meshes.size() == 1) {
            prim.material = 0;
        }

        // Indices
        if (!mesh.indices.empty()) {
            int v = ctx.appendView(mesh.indices.data(),
                                   mesh.indices.size() * sizeof(uint32_t),
                                   TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER);
            int a = ctx.addAccessor(v, TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT,
                                    TINYGLTF_TYPE_SCALAR, mesh.indices.size());
            prim.indices = a;
        }

        gltfMesh.primitives.push_back(prim);
        ctx.model.meshes.push_back(gltfMesh);
        int gltfMeshIdx = static_cast<int>(ctx.model.meshes.size()) - 1;

        tinygltf::Node meshNode;
        meshNode.mesh = gltfMeshIdx;
        if (skinned) {
            meshNode.skin = skeletonSkinIndex[skIdx];
        }
        ctx.model.nodes.push_back(meshNode);
        meshNodeIndices.push_back(static_cast<int>(ctx.model.nodes.size()) - 1);
    }

    // --- Scene root nodes ---
    tinygltf::Scene gltfScene;
    for (int nodeIdx : meshNodeIndices) {
        gltfScene.nodes.push_back(nodeIdx);
    }
    for (size_t s = 0; s < scene.skeletons.size(); ++s) {
        const auto& sk = scene.skeletons[s];
        int fbn = skeletonFirstBoneNode[s];
        for (size_t b = 0; b < sk.bones.size(); ++b) {
            if (sk.bones[b].parent < 0) {
                gltfScene.nodes.push_back(fbn + static_cast<int>(b));
            }
        }
    }
    ctx.model.scenes.push_back(gltfScene);
    ctx.model.defaultScene = 0;

    // --- Animations ---
    for (size_t ai = 0; ai < scene.animations.size(); ++ai) {
        const auto& anim = scene.animations[ai];
        int skelIdx = (ai < scene.animationSkeleton.size()) ? scene.animationSkeleton[ai] : 0;
        if (skelIdx < 0 || skelIdx >= static_cast<int>(scene.skeletons.size())) continue;
        int firstBoneNode = skeletonFirstBoneNode[skelIdx];
        const auto& sk = scene.skeletons[skelIdx];

        tinygltf::Animation gAnim;
        gAnim.name = anim.name;
        for (const auto& ch : anim.channels) {
            if (ch.boneIndex < 0 || ch.boneIndex >= static_cast<int>(sk.bones.size())) continue;

            int vIn = ctx.appendView(ch.times.data(),
                                     ch.times.size() * sizeof(float));
            float tmin = ch.times.empty() ? 0.0f : ch.times.front();
            float tmax = ch.times.empty() ? 0.0f : ch.times.back();
            int aIn = ctx.addAccessor(vIn, TINYGLTF_COMPONENT_TYPE_FLOAT,
                                      TINYGLTF_TYPE_SCALAR, ch.times.size(),
                                      {(double)tmin}, {(double)tmax});

            int vOut = ctx.appendView(ch.values.data(),
                                      ch.values.size() * sizeof(float));
            int outType = (ch.path == AnimChannel::Path::Rotation)
                              ? TINYGLTF_TYPE_VEC4 : TINYGLTF_TYPE_VEC3;
            int stride = (ch.path == AnimChannel::Path::Rotation) ? 4 : 3;
            // Count is in elements; a CUBICSPLINE sampler's output holds 3
            // elements per keyframe and the glTF count includes all three.
            size_t outCount = ch.values.size() / stride;
            int aOut = ctx.addAccessor(vOut, TINYGLTF_COMPONENT_TYPE_FLOAT,
                                       outType, outCount);

            tinygltf::AnimationSampler s;
            s.input = aIn;
            s.output = aOut;
            switch (ch.interp) {
                case AnimChannel::Interp::Step: s.interpolation = "STEP"; break;
                case AnimChannel::Interp::CubicSpline: s.interpolation = "CUBICSPLINE"; break;
                default: s.interpolation = "LINEAR"; break;
            }
            gAnim.samplers.push_back(s);

            tinygltf::AnimationChannel gCh;
            gCh.sampler = static_cast<int>(gAnim.samplers.size()) - 1;
            gCh.target_node = firstBoneNode + ch.boneIndex;
            switch (ch.path) {
                case AnimChannel::Path::Translation: gCh.target_path = "translation"; break;
                case AnimChannel::Path::Rotation:    gCh.target_path = "rotation"; break;
                case AnimChannel::Path::Scale:       gCh.target_path = "scale"; break;
            }
            gAnim.channels.push_back(gCh);
        }
        if (!gAnim.channels.empty())
            ctx.model.animations.push_back(gAnim);
    }

    if (!ctx.buffer.data.empty()) {
        ctx.model.buffers.push_back(std::move(ctx.buffer));
    }

    tinygltf::TinyGLTF writer;
    bool binary = path.size() >= 4 && path.substr(path.size() - 4) == ".glb";
    return writer.WriteGltfSceneToFile(&ctx.model, path, true, true, true, binary);
}

bool saveGLTF(const MeshData& mesh, const std::string& path) {
    if (mesh.empty()) return false;
    GltfScene scene;
    scene.meshes.push_back(mesh);
    scene.meshSkeleton.push_back(-1);
    scene.meshMaterial.push_back(-1);
    return saveGLTF(scene, path);
}

bool saveGLTF(const MeshData& mesh,
              const SkinData* skin,
              const Skeleton* skeleton,
              const std::vector<Animation>& animations,
              const std::string& path) {
    if (mesh.empty()) return false;
    GltfScene scene;
    scene.meshes.push_back(mesh);
    if (skin) scene.skins.push_back(*skin);
    else scene.skins.push_back({});
    if (skeleton) {
        scene.skeletons.push_back(*skeleton);
        scene.meshSkeleton.push_back(0);
    } else {
        scene.meshSkeleton.push_back(-1);
    }
    scene.meshMaterial.push_back(-1);
    scene.animations = animations;
    scene.animationSkeleton.resize(animations.size(), 0);
    return saveGLTF(scene, path);
}

bool saveGLTF(const MeshData& mesh,
              std::optional<std::reference_wrapper<const SkinData>> skin,
              std::optional<std::reference_wrapper<const Skeleton>> skeleton,
              const std::vector<Animation>& animations,
              const std::string& path) {
    const SkinData* sk = skin ? &skin->get() : nullptr;
    const Skeleton* skel = skeleton ? &skeleton->get() : nullptr;
    return saveGLTF(mesh, sk, skel, animations, path);
}

#endif // BROMESH_HAS_GLTF

} // namespace bromesh
