#include "gltf_common.h"

#if BROMESH_HAS_DRACO
#include "bromesh/io/draco.h"
#endif

namespace bromesh {

#if !BROMESH_HAS_GLTF
GltfScene loadGLTF(const std::string&) { return {}; }
#else

#if BROMESH_HAS_DRACO
static bool extractDracoFloats(const DracoAttribute& attr, int targetComps, std::vector<float>& out) {
    if (attr.count == 0 || attr.components <= 0) return false;
    out.resize(static_cast<size_t>(attr.count) * targetComps, 0.0f);
    int comps = std::min(attr.components, targetComps);
    for (uint32_t i = 0; i < attr.count; ++i) {
        for (int c = 0; c < comps; ++c) {
            size_t idx = static_cast<size_t>(i) * attr.components + c;
            float val = 0.0f;
            switch (attr.kind) {
                case DracoAttribute::Kind::Float32:
                    val = reinterpret_cast<const float*>(attr.bytes.data())[idx];
                    break;
                case DracoAttribute::Kind::Uint8:
                    val = reinterpret_cast<const uint8_t*>(attr.bytes.data())[idx] / 255.0f;
                    break;
                case DracoAttribute::Kind::Uint16:
                    val = reinterpret_cast<const uint16_t*>(attr.bytes.data())[idx] / 65535.0f;
                    break;
                case DracoAttribute::Kind::Int8:
                    val = std::max(-1.0f, reinterpret_cast<const int8_t*>(attr.bytes.data())[idx] / 127.0f);
                    break;
                case DracoAttribute::Kind::Int16:
                    val = std::max(-1.0f, reinterpret_cast<const int16_t*>(attr.bytes.data())[idx] / 32767.0f);
                    break;
                default: break;
            }
            out[static_cast<size_t>(i) * targetComps + c] = val;
        }
        if (targetComps == 4 && comps == 3) {
            out[static_cast<size_t>(i) * 4 + 3] = 1.0f;
        }
    }
    return true;
}

static bool extractDracoUints(const DracoAttribute& attr, int targetComps, std::vector<uint32_t>& out) {
    if (attr.count == 0 || attr.components <= 0) return false;
    out.resize(static_cast<size_t>(attr.count) * targetComps, 0);
    int comps = std::min(attr.components, targetComps);
    for (uint32_t i = 0; i < attr.count; ++i) {
        for (int c = 0; c < comps; ++c) {
            size_t idx = static_cast<size_t>(i) * attr.components + c;
            uint32_t val = 0;
            switch (attr.kind) {
                case DracoAttribute::Kind::Uint8:
                    val = reinterpret_cast<const uint8_t*>(attr.bytes.data())[idx];
                    break;
                case DracoAttribute::Kind::Uint16:
                    val = reinterpret_cast<const uint16_t*>(attr.bytes.data())[idx];
                    break;
                case DracoAttribute::Kind::Uint32:
                    val = reinterpret_cast<const uint32_t*>(attr.bytes.data())[idx];
                    break;
                case DracoAttribute::Kind::Int8:
                    val = static_cast<uint32_t>(std::max(0, static_cast<int>(reinterpret_cast<const int8_t*>(attr.bytes.data())[idx])));
                    break;
                case DracoAttribute::Kind::Int16:
                    val = static_cast<uint32_t>(std::max(0, static_cast<int>(reinterpret_cast<const int16_t*>(attr.bytes.data())[idx])));
                    break;
                case DracoAttribute::Kind::Int32:
                    val = static_cast<uint32_t>(std::max(0, reinterpret_cast<const int32_t*>(attr.bytes.data())[idx]));
                    break;
                default: break;
            }
            out[static_cast<size_t>(i) * targetComps + c] = val;
        }
    }
    return true;
}
#endif

GltfScene loadGLTF(const std::string& path) {
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

    // Per-skin: map of glTF node index -> bone index within that skeleton.
    std::vector<std::unordered_map<int, int>> nodeToBone(model.skins.size());

    // Build a parent-of map for all nodes
    std::vector<int> parentOfNode(model.nodes.size(), -1);
    for (size_t n = 0; n < model.nodes.size(); ++n) {
        for (int c : model.nodes[n].children) {
            if (c >= 0 && static_cast<size_t>(c) < parentOfNode.size()) {
                parentOfNode[c] = static_cast<int>(n);
            }
        }
    }

    // Build skeletons from glTF skins.
    for (size_t si = 0; si < model.skins.size(); ++si) {
        const auto& gskin = model.skins[si];
        Skeleton sk;
        sk.bones.resize(gskin.joints.size());

        auto& n2b = nodeToBone[si];
        for (size_t j = 0; j < gskin.joints.size(); ++j)
            n2b[gskin.joints[j]] = static_cast<int>(j);

        std::vector<float> ibm;
        if (gskin.inverseBindMatrices >= 0) {
            const auto& acc = model.accessors[gskin.inverseBindMatrices];
            const auto* src = reinterpret_cast<const float*>(accessorData(model, gskin.inverseBindMatrices));
            if (src) ibm.assign(src, src + acc.count * 16);
        }

        for (size_t j = 0; j < gskin.joints.size(); ++j) {
            int nodeIdx = gskin.joints[j];
            const auto& node = model.nodes[nodeIdx];
            Bone& b = sk.bones[j];
            b.name = node.name;

            b.parent = -1;
            for (size_t n = 0; n < model.nodes.size(); ++n) {
                const auto& cand = model.nodes[n];
                for (int c : cand.children) {
                    if (c == nodeIdx) {
                        auto it = n2b.find(static_cast<int>(n));
                        if (it != n2b.end()) b.parent = it->second;
                        break;
                    }
                }
                if (b.parent != -1) break;
            }

            if (!node.matrix.empty()) {
                float m[16];
                for (int k = 0; k < 16; ++k) m[k] = static_cast<float>(node.matrix[k]);
                decomposeTRS(m, b.localT, b.localR, b.localS);
            } else {
                if (node.translation.size() == 3) {
                    b.localT[0] = static_cast<float>(node.translation[0]);
                    b.localT[1] = static_cast<float>(node.translation[1]);
                    b.localT[2] = static_cast<float>(node.translation[2]);
                }
                if (node.rotation.size() == 4) {
                    b.localR[0] = static_cast<float>(node.rotation[0]);
                    b.localR[1] = static_cast<float>(node.rotation[1]);
                    b.localR[2] = static_cast<float>(node.rotation[2]);
                    b.localR[3] = static_cast<float>(node.rotation[3]);
                }
                if (node.scale.size() == 3) {
                    b.localS[0] = static_cast<float>(node.scale[0]);
                    b.localS[1] = static_cast<float>(node.scale[1]);
                    b.localS[2] = static_cast<float>(node.scale[2]);
                }
            }

            if (!ibm.empty() && j * 16 + 16 <= ibm.size()) {
                std::memcpy(b.inverseBind, &ibm[j * 16], 16 * sizeof(float));
            } else {
                matIdentity(b.inverseBind);
            }
        }

        for (size_t j = 0; j < sk.bones.size(); ++j) {
            if (sk.bones[j].parent != -1) continue;

            int nodeIdx = gskin.joints[j];
            float anc[16]; matIdentity(anc);
            std::vector<int> chain;
            for (int cur = parentOfNode[nodeIdx]; cur != -1; cur = parentOfNode[cur]) {
                chain.push_back(cur);
            }
            for (auto it = chain.rbegin(); it != chain.rend(); ++it) {
                float lm[16]; nodeLocalMat(model.nodes[*it], lm);
                float out[16]; matMul4(anc, lm, out);
                std::memcpy(anc, out, 16 * sizeof(float));
            }
            std::memcpy(sk.rootTransform, anc, 16 * sizeof(float));
            break;
        }

        scene.skeletons.push_back(std::move(sk));
    }

    std::vector<int> meshSkinIndex(model.meshes.size(), -1);
    for (const auto& node : model.nodes) {
        if (node.mesh >= 0 && node.skin >= 0 && static_cast<size_t>(node.mesh) < meshSkinIndex.size())
            meshSkinIndex[node.mesh] = node.skin;
    }

    // Build meshes
    for (size_t mi = 0; mi < model.meshes.size(); ++mi) {
        const auto& mesh = model.meshes[mi];
        for (const auto& prim : mesh.primitives) {
            if (prim.mode != TINYGLTF_MODE_TRIANGLES) continue;

            MeshData md;
            SkinData skin;
            bool loadedViaDraco = false;

#if BROMESH_HAS_DRACO
            auto dracoIt = prim.extensions.find("KHR_draco_mesh_compression");
            if (dracoIt != prim.extensions.end()) {
                const auto& ext = dracoIt->second;
                if (ext.IsObject() && ext.Has("bufferView")) {
                    int bvIdx = ext.Get("bufferView").GetNumberAsInt();
                    if (bvIdx >= 0 && bvIdx < static_cast<int>(model.bufferViews.size())) {
                        const auto& bv = model.bufferViews[bvIdx];
                        if (bv.buffer >= 0 && bv.buffer < static_cast<int>(model.buffers.size())) {
                            const auto& buf = model.buffers[bv.buffer];
                            if (bv.byteOffset + bv.byteLength <= buf.data.size()) {
                                const uint8_t* compressedData = buf.data.data() + bv.byteOffset;
                                DracoDecoded decoded = decodeDraco(compressedData, bv.byteLength);
                                if (decoded.ok()) {
                                    loadedViaDraco = true;
                                    auto findAttr = [&](uint32_t id) -> const DracoAttribute* {
                                        for (const auto& a : decoded.attributes) {
                                            if (a.uniqueId == id) return &a;
                                        }
                                        return nullptr;
                                    };

                                    const tinygltf::Value& attrs = ext.Has("attributes") ? ext.Get("attributes") : tinygltf::Value();

                                    // Positions
                                    if (attrs.IsObject() && attrs.Has("POSITION")) {
                                        if (!decoded.mesh.positions.empty()) {
                                            md.positions = std::move(decoded.mesh.positions);
                                        } else {
                                            uint32_t id = static_cast<uint32_t>(attrs.Get("POSITION").GetNumberAsInt());
                                            if (const auto* a = findAttr(id)) {
                                                extractDracoFloats(*a, 3, md.positions);
                                            }
                                        }
                                    } else if (!decoded.mesh.positions.empty()) {
                                        md.positions = std::move(decoded.mesh.positions);
                                    }

                                    // Indices
                                    md.indices = std::move(decoded.mesh.indices);
                                    if (md.indices.empty() && !md.positions.empty()) {
                                        size_t n = md.positions.size() / 3;
                                        md.indices.resize(n);
                                        for (size_t i = 0; i < n; ++i) md.indices[i] = static_cast<uint32_t>(i);
                                    }

                                    // Normals
                                    if (attrs.IsObject() && attrs.Has("NORMAL")) {
                                        if (!decoded.mesh.normals.empty()) {
                                            md.normals = std::move(decoded.mesh.normals);
                                        } else {
                                            uint32_t id = static_cast<uint32_t>(attrs.Get("NORMAL").GetNumberAsInt());
                                            if (const auto* a = findAttr(id)) {
                                                extractDracoFloats(*a, 3, md.normals);
                                            }
                                        }
                                    } else if (!decoded.mesh.normals.empty()) {
                                        md.normals = std::move(decoded.mesh.normals);
                                    }

                                    // UVs
                                    if (attrs.IsObject() && attrs.Has("TEXCOORD_0")) {
                                        if (!decoded.mesh.uvs.empty()) {
                                            md.uvs = std::move(decoded.mesh.uvs);
                                        } else {
                                            uint32_t id = static_cast<uint32_t>(attrs.Get("TEXCOORD_0").GetNumberAsInt());
                                            if (const auto* a = findAttr(id)) {
                                                extractDracoFloats(*a, 2, md.uvs);
                                            }
                                        }
                                    } else if (!decoded.mesh.uvs.empty()) {
                                        md.uvs = std::move(decoded.mesh.uvs);
                                    }

                                    // Tangents
                                    if (attrs.IsObject() && attrs.Has("TANGENT")) {
                                        uint32_t id = static_cast<uint32_t>(attrs.Get("TANGENT").GetNumberAsInt());
                                        if (const auto* a = findAttr(id)) {
                                            extractDracoFloats(*a, 4, md.tangents);
                                        }
                                    }

                                    // Colors
                                    if (attrs.IsObject() && attrs.Has("COLOR_0")) {
                                        if (!decoded.mesh.colors.empty()) {
                                            md.colors = std::move(decoded.mesh.colors);
                                        } else {
                                            uint32_t id = static_cast<uint32_t>(attrs.Get("COLOR_0").GetNumberAsInt());
                                            if (const auto* a = findAttr(id)) {
                                                extractDracoFloats(*a, 4, md.colors);
                                            }
                                        }
                                    } else if (!decoded.mesh.colors.empty()) {
                                        md.colors = std::move(decoded.mesh.colors);
                                    }

                                    // Skinning: JOINTS_0 and WEIGHTS_0
                                    if (attrs.IsObject() && attrs.Has("JOINTS_0") && attrs.Has("WEIGHTS_0")) {
                                        uint32_t jId = static_cast<uint32_t>(attrs.Get("JOINTS_0").GetNumberAsInt());
                                        uint32_t wId = static_cast<uint32_t>(attrs.Get("WEIGHTS_0").GetNumberAsInt());
                                        const auto* jAttr = findAttr(jId);
                                        const auto* wAttr = findAttr(wId);
                                        if (jAttr && wAttr) {
                                            extractDracoUints(*jAttr, 4, skin.boneIndices);
                                            extractDracoFloats(*wAttr, 4, skin.boneWeights);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
#endif // BROMESH_HAS_DRACO

            if (!loadedViaDraco) {
                // Indices
                if (prim.indices >= 0) {
                    const auto& accessor = model.accessors[prim.indices];
                    const uint8_t* base = accessorData(model, prim.indices);
                    if (base) {
                        md.indices.resize(accessor.count);
                        if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT) {
                            const auto* src = reinterpret_cast<const uint16_t*>(base);
                            for (size_t i = 0; i < accessor.count; ++i) md.indices[i] = src[i];
                        } else if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT) {
                            const auto* src = reinterpret_cast<const uint32_t*>(base);
                            for (size_t i = 0; i < accessor.count; ++i) md.indices[i] = src[i];
                        } else if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE) {
                            const uint8_t* src = base;
                            for (size_t i = 0; i < accessor.count; ++i) md.indices[i] = src[i];
                        }
                    }
                }

                auto readVec = [&](const char* attrName, std::vector<float>& dst, int components) {
                    auto it = prim.attributes.find(attrName);
                    if (it == prim.attributes.end()) return;
                    const auto& accessor = model.accessors[it->second];
                    const auto* src = reinterpret_cast<const float*>(accessorData(model, it->second));
                    if (src) dst.assign(src, src + accessor.count * components);
                };

                readVec("POSITION", md.positions, 3);
                readVec("NORMAL",   md.normals,   3);
                readVec("TEXCOORD_0", md.uvs,     2);
                readVec("TANGENT",  md.tangents,  4);

                if (prim.indices < 0 && !md.positions.empty()) {
                    size_t n = md.positions.size() / 3;
                    md.indices.resize(n);
                    for (size_t i = 0; i < n; ++i) md.indices[i] = static_cast<uint32_t>(i);
                }

                // COLOR_0
                {
                    auto it = prim.attributes.find("COLOR_0");
                    if (it != prim.attributes.end()) {
                        const auto& accessor = model.accessors[it->second];
                        if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT) {
                            int comps = (accessor.type == TINYGLTF_TYPE_VEC4) ? 4 : 3;
                            const auto* src = reinterpret_cast<const float*>(accessorData(model, it->second));
                            if (src) {
                                md.colors.resize(accessor.count * 4);
                                for (size_t i = 0; i < accessor.count; ++i) {
                                    md.colors[i*4 + 0] = src[i*comps + 0];
                                    md.colors[i*4 + 1] = src[i*comps + 1];
                                    md.colors[i*4 + 2] = src[i*comps + 2];
                                    md.colors[i*4 + 3] = (comps == 4) ? src[i*comps + 3] : 1.0f;
                                }
                            }
                        }
                    }
                }

                // Joints/weights
                {
                    auto itJ = prim.attributes.find("JOINTS_0");
                    auto itW = prim.attributes.find("WEIGHTS_0");
                    if (itJ != prim.attributes.end() && itW != prim.attributes.end()) {
                        const auto& accJ = model.accessors[itJ->second];
                        const auto& accW = model.accessors[itW->second];
                        skin.boneIndices.resize(accJ.count * 4);
                        skin.boneWeights.resize(accW.count * 4);

                        const uint8_t* jbase = accessorData(model, itJ->second);
                        if (jbase) {
                            if (accJ.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT) {
                                const auto* src = reinterpret_cast<const uint16_t*>(jbase);
                                for (size_t i = 0; i < accJ.count * 4; ++i) skin.boneIndices[i] = src[i];
                            } else {
                                for (size_t i = 0; i < accJ.count * 4; ++i) skin.boneIndices[i] = jbase[i];
                            }
                        }

                        if (accW.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT) {
                            const auto* src = reinterpret_cast<const float*>(accessorData(model, itW->second));
                            if (src) {
                                std::memcpy(skin.boneWeights.data(), src, accW.count * 4 * sizeof(float));
                            }
                        }
                    }
                }
            }

            int skIdx = meshSkinIndex[mi];
            if (skIdx >= 0 && skIdx < static_cast<int>(scene.skeletons.size())) {
                skin.boneCount = scene.skeletons[skIdx].bones.size();
                skin.inverseBindMatrices.resize(skin.boneCount * 16);
                for (size_t j = 0; j < skin.boneCount; ++j)
                    std::memcpy(&skin.inverseBindMatrices[j * 16],
                                scene.skeletons[skIdx].bones[j].inverseBind,
                                16 * sizeof(float));
            }

            scene.meshes.push_back(std::move(md));
            scene.skins.push_back(std::move(skin));
            scene.meshSkeleton.push_back(skIdx);
            scene.meshMaterial.push_back(prim.material);
        }
    }

    // Build images
    scene.images.reserve(model.images.size());
    for (const auto& gi : model.images) {
        Image img;
        img.name = gi.name;
        img.mimeType = gi.mimeType;
        img.width  = gi.width;
        img.height = gi.height;

        if (gi.width > 0 && gi.height > 0 && !gi.image.empty()) {
            const size_t px = static_cast<size_t>(gi.width) * static_cast<size_t>(gi.height);
            img.data.resize(px * 4);
            const int c = gi.component;
            const int bits = gi.bits;
            const uint8_t* src = gi.image.data();

            auto get8 = [&](size_t i, int channel) -> uint8_t {
                if (bits == 16) {
                    uint16_t v = (reinterpret_cast<const uint16_t*>(src))[i * c + channel];
                    return static_cast<uint8_t>(v >> 8);
                }
                return src[i * c + channel];
            };
            for (size_t i = 0; i < px; ++i) {
                uint8_t r = 0, g = 0, b = 0, a = 255;
                if (c >= 1) r = get8(i, 0);
                if (c >= 2) g = get8(i, 1); else g = r;
                if (c >= 3) b = get8(i, 2); else b = r;
                if (c >= 4) a = get8(i, 3);
                img.data[i*4 + 0] = r;
                img.data[i*4 + 1] = g;
                img.data[i*4 + 2] = b;
                img.data[i*4 + 3] = a;
            }
        }
        scene.images.push_back(std::move(img));
    }

    // Build materials
    scene.materials.reserve(model.materials.size());
    auto resolveImage = [&](int texIdx) -> int {
        if (texIdx < 0 || texIdx >= static_cast<int>(model.textures.size())) return -1;
        int src = model.textures[texIdx].source;
        if (src < 0 || src >= static_cast<int>(scene.images.size())) return -1;
        return src;
    };
    for (const auto& gm : model.materials) {
        Material mat;
        mat.name = gm.name;
        const auto& pbr = gm.pbrMetallicRoughness;
        if (pbr.baseColorFactor.size() == 4) {
            mat.baseColorFactor[0] = static_cast<float>(pbr.baseColorFactor[0]);
            mat.baseColorFactor[1] = static_cast<float>(pbr.baseColorFactor[1]);
            mat.baseColorFactor[2] = static_cast<float>(pbr.baseColorFactor[2]);
            mat.baseColorFactor[3] = static_cast<float>(pbr.baseColorFactor[3]);
        }
        mat.metallicFactor  = static_cast<float>(pbr.metallicFactor);
        mat.roughnessFactor = static_cast<float>(pbr.roughnessFactor);
        if (gm.emissiveFactor.size() == 3) {
            mat.emissiveFactor[0] = static_cast<float>(gm.emissiveFactor[0]);
            mat.emissiveFactor[1] = static_cast<float>(gm.emissiveFactor[1]);
            mat.emissiveFactor[2] = static_cast<float>(gm.emissiveFactor[2]);
        }
        mat.baseColorTexture         = resolveImage(pbr.baseColorTexture.index);
        mat.metallicRoughnessTexture = resolveImage(pbr.metallicRoughnessTexture.index);
        mat.normalTexture            = resolveImage(gm.normalTexture.index);
        mat.occlusionTexture         = resolveImage(gm.occlusionTexture.index);
        mat.emissiveTexture          = resolveImage(gm.emissiveTexture.index);
        scene.materials.push_back(std::move(mat));
    }

    // Build animations
    for (const auto& ganim : model.animations) {
        Animation anim;
        anim.name = ganim.name;
        anim.duration = 0.0f;
        int associatedSkel = -1;

        for (const auto& gch : ganim.channels) {
            AnimChannel ch;

            int boneIdx = -1;
            int skelIdx = -1;
            for (size_t si = 0; si < nodeToBone.size(); ++si) {
                auto it = nodeToBone[si].find(gch.target_node);
                if (it != nodeToBone[si].end()) {
                    boneIdx = it->second;
                    skelIdx = static_cast<int>(si);
                    break;
                }
            }
            if (boneIdx < 0) continue;
            if (associatedSkel < 0) associatedSkel = skelIdx;

            ch.boneIndex = boneIdx;
            if (gch.target_path == "translation") ch.path = AnimChannel::Path::Translation;
            else if (gch.target_path == "rotation") ch.path = AnimChannel::Path::Rotation;
            else if (gch.target_path == "scale")    ch.path = AnimChannel::Path::Scale;
            else continue;

            const auto& gs = ganim.samplers[gch.sampler];
            if (gs.interpolation == "STEP") ch.interp = AnimChannel::Interp::Step;
            else if (gs.interpolation == "CUBICSPLINE") ch.interp = AnimChannel::Interp::CubicSpline;
            else ch.interp = AnimChannel::Interp::Linear;

            const auto& accIn = model.accessors[gs.input];
            const auto* tsrc = reinterpret_cast<const float*>(accessorData(model, gs.input));
            if (tsrc) {
                ch.times.assign(tsrc, tsrc + accIn.count);
                if (!ch.times.empty()) {
                    anim.duration = std::max(anim.duration, ch.times.back());
                }
            }

            const auto& accOut = model.accessors[gs.output];
            const auto* vsrc = reinterpret_cast<const float*>(accessorData(model, gs.output));
            int stride = (ch.path == AnimChannel::Path::Rotation) ? 4 : 3;
            int packing = (ch.interp == AnimChannel::Interp::CubicSpline) ? 3 : 1;
            if (vsrc) {
                ch.values.assign(vsrc, vsrc + accOut.count * stride * packing);
            }

            anim.channels.push_back(std::move(ch));
        }

        scene.animations.push_back(std::move(anim));
        scene.animationSkeleton.push_back(associatedSkel);
    }

    return scene;
}

#endif // BROMESH_HAS_GLTF

} // namespace bromesh
