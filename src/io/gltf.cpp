#include "bromesh/io/gltf.h"

#if BROMESH_HAS_GLTF
#include "tiny_gltf.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <unordered_map>
#endif

namespace bromesh {

#if BROMESH_HAS_GLTF

// ---- small helpers ---------------------------------------------------------

static void matIdentity(float* m) {
    for (int i = 0; i < 16; ++i) m[i] = 0.0f;
    m[0] = m[5] = m[10] = m[15] = 1.0f;
}

// Decompose a column-major 4x4 affine matrix into T (vec3), R (quat xyzw),
// S (vec3). Assumes affine with no shear.
static void decomposeTRS(const float* m, float t[3], float r[4], float s[3]) {
    t[0] = m[12]; t[1] = m[13]; t[2] = m[14];

    float cx[3] = { m[0], m[1], m[2] };
    float cy[3] = { m[4], m[5], m[6] };
    float cz[3] = { m[8], m[9], m[10] };

    s[0] = std::sqrt(cx[0]*cx[0] + cx[1]*cx[1] + cx[2]*cx[2]);
    s[1] = std::sqrt(cy[0]*cy[0] + cy[1]*cy[1] + cy[2]*cy[2]);
    s[2] = std::sqrt(cz[0]*cz[0] + cz[1]*cz[1] + cz[2]*cz[2]);

    if (s[0] == 0 || s[1] == 0 || s[2] == 0) {
        r[0] = 0; r[1] = 0; r[2] = 0; r[3] = 1;
        return;
    }

    float rot[9] = {
        cx[0]/s[0], cx[1]/s[0], cx[2]/s[0],
        cy[0]/s[1], cy[1]/s[1], cy[2]/s[1],
        cz[0]/s[2], cz[1]/s[2], cz[2]/s[2],
    };
    // rot is column-major 3x3. Convert to quaternion (xyzw).
    float trace = rot[0] + rot[4] + rot[8];
    if (trace > 0.0f) {
        float k = std::sqrt(trace + 1.0f) * 2.0f;
        r[3] = 0.25f * k;
        r[0] = (rot[5] - rot[7]) / k;
        r[1] = (rot[6] - rot[2]) / k;
        r[2] = (rot[1] - rot[3]) / k;
    } else if (rot[0] > rot[4] && rot[0] > rot[8]) {
        float k = std::sqrt(1.0f + rot[0] - rot[4] - rot[8]) * 2.0f;
        r[3] = (rot[5] - rot[7]) / k;
        r[0] = 0.25f * k;
        r[1] = (rot[3] + rot[1]) / k;
        r[2] = (rot[6] + rot[2]) / k;
    } else if (rot[4] > rot[8]) {
        float k = std::sqrt(1.0f + rot[4] - rot[0] - rot[8]) * 2.0f;
        r[3] = (rot[6] - rot[2]) / k;
        r[0] = (rot[3] + rot[1]) / k;
        r[1] = 0.25f * k;
        r[2] = (rot[7] + rot[5]) / k;
    } else {
        float k = std::sqrt(1.0f + rot[8] - rot[0] - rot[4]) * 2.0f;
        r[3] = (rot[1] - rot[3]) / k;
        r[0] = (rot[6] + rot[2]) / k;
        r[1] = (rot[7] + rot[5]) / k;
        r[2] = 0.25f * k;
    }
}

// Pointer to raw accessor data.
static const uint8_t* accessorData(const tinygltf::Model& model, int accessorIdx) {
    const auto& a = model.accessors[accessorIdx];
    const auto& bv = model.bufferViews[a.bufferView];
    return model.buffers[bv.buffer].data.data() + bv.byteOffset + a.byteOffset;
}

// Column-major 4x4 multiply: out = a * b.
static void matMul4(const float* a, const float* b, float* out) {
    float tmp[16];
    for (int c = 0; c < 4; ++c) {
        for (int r = 0; r < 4; ++r) {
            tmp[c*4 + r] =
                a[0*4+r]*b[c*4+0] + a[1*4+r]*b[c*4+1] +
                a[2*4+r]*b[c*4+2] + a[3*4+r]*b[c*4+3];
        }
    }
    std::memcpy(out, tmp, 16 * sizeof(float));
}

// Build a column-major 4x4 from T, R (quat xyzw), S.
static void trsToMat(const float t[3], const float r[4], const float s[3], float* m) {
    float xx = r[0]*r[0], yy = r[1]*r[1], zz = r[2]*r[2];
    float xy = r[0]*r[1], xz = r[0]*r[2], yz = r[1]*r[2];
    float wx = r[3]*r[0], wy = r[3]*r[1], wz = r[3]*r[2];
    m[0]  = (1.0f - 2.0f*(yy + zz)) * s[0];
    m[1]  = (2.0f*(xy + wz))        * s[0];
    m[2]  = (2.0f*(xz - wy))        * s[0];
    m[3]  = 0.0f;
    m[4]  = (2.0f*(xy - wz))        * s[1];
    m[5]  = (1.0f - 2.0f*(xx + zz)) * s[1];
    m[6]  = (2.0f*(yz + wx))        * s[1];
    m[7]  = 0.0f;
    m[8]  = (2.0f*(xz + wy))        * s[2];
    m[9]  = (2.0f*(yz - wx))        * s[2];
    m[10] = (1.0f - 2.0f*(xx + yy)) * s[2];
    m[11] = 0.0f;
    m[12] = t[0]; m[13] = t[1]; m[14] = t[2]; m[15] = 1.0f;
}

// Build the local transform matrix of a glTF node (matrix, or TRS).
static void nodeLocalMat(const tinygltf::Node& node, float* m) {
    if (!node.matrix.empty() && node.matrix.size() == 16) {
        for (int k = 0; k < 16; ++k) m[k] = static_cast<float>(node.matrix[k]);
        return;
    }
    float t[3] = {0,0,0}, r[4] = {0,0,0,1}, s[3] = {1,1,1};
    if (node.translation.size() == 3) {
        t[0] = (float)node.translation[0];
        t[1] = (float)node.translation[1];
        t[2] = (float)node.translation[2];
    }
    if (node.rotation.size() == 4) {
        r[0] = (float)node.rotation[0];
        r[1] = (float)node.rotation[1];
        r[2] = (float)node.rotation[2];
        r[3] = (float)node.rotation[3];
    }
    if (node.scale.size() == 3) {
        s[0] = (float)node.scale[0];
        s[1] = (float)node.scale[1];
        s[2] = (float)node.scale[2];
    }
    trsToMat(t, r, s, m);
}

// ---- load ------------------------------------------------------------------

GltfScene loadGLTF(const std::string& path) {
    tinygltf::Model model;
    tinygltf::TinyGLTF loader;
    std::string err, warn;

    // Let tinygltf's bundled stb_image decode embedded textures.
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

    // Build a parent-of map for all nodes (node idx -> parent node idx, -1 if
    // none). glTF nodes only know their children, so we invert the relation.
    std::vector<int> parentOfNode(model.nodes.size(), -1);
    for (size_t n = 0; n < model.nodes.size(); ++n) {
        for (int c : model.nodes[n].children) {
            if (c >= 0 && (size_t)c < parentOfNode.size()) {
                parentOfNode[c] = static_cast<int>(n);
            }
        }
    }

    // Build skeletons from glTF skins.
    for (size_t si = 0; si < model.skins.size(); ++si) {
        const auto& gskin = model.skins[si];
        Skeleton sk;
        sk.bones.resize(gskin.joints.size());

        // node index -> local joint index
        auto& n2b = nodeToBone[si];
        for (size_t j = 0; j < gskin.joints.size(); ++j)
            n2b[gskin.joints[j]] = static_cast<int>(j);

        // Inverse bind matrices
        std::vector<float> ibm;
        if (gskin.inverseBindMatrices >= 0) {
            const auto& acc = model.accessors[gskin.inverseBindMatrices];
            const auto* src = reinterpret_cast<const float*>(accessorData(model, gskin.inverseBindMatrices));
            ibm.assign(src, src + acc.count * 16);
        }

        for (size_t j = 0; j < gskin.joints.size(); ++j) {
            int nodeIdx = gskin.joints[j];
            const auto& node = model.nodes[nodeIdx];
            Bone& b = sk.bones[j];
            b.name = node.name;

            // Parent: walk all nodes' children lists; whichever is parent and
            // also in the joints list becomes the parent bone.
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

            // Local TRS
            if (!node.matrix.empty()) {
                float m[16];
                for (int k = 0; k < 16; ++k) m[k] = static_cast<float>(node.matrix[k]);
                decomposeTRS(m, b.localT, b.localR, b.localS);
            } else {
                if (node.translation.size() == 3) {
                    b.localT[0] = (float)node.translation[0];
                    b.localT[1] = (float)node.translation[1];
                    b.localT[2] = (float)node.translation[2];
                }
                if (node.rotation.size() == 4) {
                    b.localR[0] = (float)node.rotation[0];
                    b.localR[1] = (float)node.rotation[1];
                    b.localR[2] = (float)node.rotation[2];
                    b.localR[3] = (float)node.rotation[3];
                }
                if (node.scale.size() == 3) {
                    b.localS[0] = (float)node.scale[0];
                    b.localS[1] = (float)node.scale[1];
                    b.localS[2] = (float)node.scale[2];
                }
            }

            if (!ibm.empty() && j * 16 + 16 <= ibm.size()) {
                std::memcpy(b.inverseBind, &ibm[j * 16], 16 * sizeof(float));
            } else {
                matIdentity(b.inverseBind);
            }
        }

        // Capture any scene-tree ancestry above the joint list into
        // Skeleton::rootTransform. glTF stores inverseBindMatrices relative
        // to the full scene-root → joint chain, and animation channels are
        // authored in the armature's local space — so the armature ancestor
        // must participate in pose→world evaluation for both mesh vertices
        // (authored in world/meters) and animation translations (authored
        // in armature-local units) to land in a consistent space.
        // `computeWorldMatrices` premultiplies rootTransform onto root bones,
        // so IBMs are left untouched from the glTF file.
        //
        // All root bones in a glTF skin typically share the same armature
        // ancestor; we capture anc from the first root and apply uniformly.
        for (size_t j = 0; j < sk.bones.size(); ++j) {
            if (sk.bones[j].parent != -1) continue;

            int nodeIdx = gskin.joints[j];
            float anc[16]; matIdentity(anc);
            std::vector<int> chain;
            for (int cur = parentOfNode[nodeIdx]; cur != -1; cur = parentOfNode[cur]) {
                chain.push_back(cur);
            }
            // Walk from outermost ancestor inward so anc accumulates
            // left-to-right (column-major): anc = outermost * ... * innermost.
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

    // Build a node -> mesh-usage map: which nodes reference a mesh, and with
    // which skin. glTF meshes can be instanced; we flatten by mesh index and
    // assume a single skin per mesh.
    std::vector<int> meshSkinIndex(model.meshes.size(), -1);
    for (const auto& node : model.nodes) {
        if (node.mesh >= 0 && node.skin >= 0)
            meshSkinIndex[node.mesh] = node.skin;
    }

    // Build meshes.
    for (size_t mi = 0; mi < model.meshes.size(); ++mi) {
        const auto& mesh = model.meshes[mi];
        for (const auto& prim : mesh.primitives) {
            if (prim.mode != TINYGLTF_MODE_TRIANGLES) continue;

            MeshData md;
            SkinData skin;

            // Indices
            if (prim.indices >= 0) {
                const auto& accessor = model.accessors[prim.indices];
                const uint8_t* base = accessorData(model, prim.indices);
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

            size_t vertCount = 0;

            auto readVec = [&](const char* attrName, std::vector<float>& dst, int components) {
                auto it = prim.attributes.find(attrName);
                if (it == prim.attributes.end()) return;
                const auto& accessor = model.accessors[it->second];
                const auto* src = reinterpret_cast<const float*>(accessorData(model, it->second));
                dst.assign(src, src + accessor.count * components);
                vertCount = accessor.count;
            };

            readVec("POSITION", md.positions, 3);
            readVec("NORMAL",   md.normals,   3);
            readVec("TEXCOORD_0", md.uvs,     2);

            // COLOR_0 can be float or unsigned short / byte; handle float-only.
            {
                auto it = prim.attributes.find("COLOR_0");
                if (it != prim.attributes.end()) {
                    const auto& accessor = model.accessors[it->second];
                    if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT) {
                        int comps = (accessor.type == TINYGLTF_TYPE_VEC4) ? 4 : 3;
                        const auto* src = reinterpret_cast<const float*>(accessorData(model, it->second));
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

            // Joints/weights
            {
                auto itJ = prim.attributes.find("JOINTS_0");
                auto itW = prim.attributes.find("WEIGHTS_0");
                if (itJ != prim.attributes.end() && itW != prim.attributes.end()) {
                    const auto& accJ = model.accessors[itJ->second];
                    const auto& accW = model.accessors[itW->second];
                    skin.boneIndices.resize(accJ.count * 4);
                    skin.boneWeights.resize(accW.count * 4);

                    // JOINTS_0: uint8 or uint16
                    const uint8_t* jbase = accessorData(model, itJ->second);
                    if (accJ.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT) {
                        const auto* src = reinterpret_cast<const uint16_t*>(jbase);
                        for (size_t i = 0; i < accJ.count * 4; ++i) skin.boneIndices[i] = src[i];
                    } else {
                        for (size_t i = 0; i < accJ.count * 4; ++i) skin.boneIndices[i] = jbase[i];
                    }

                    // WEIGHTS_0: float (most common)
                    if (accW.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT) {
                        const auto* src = reinterpret_cast<const float*>(accessorData(model, itW->second));
                        std::memcpy(skin.boneWeights.data(), src, accW.count * 4 * sizeof(float));
                    }
                }
            }

            // Associate skeleton
            int skIdx = meshSkinIndex[mi];
            if (skIdx >= 0 && skIdx < (int)scene.skeletons.size()) {
                skin.boneCount = scene.skeletons[skIdx].bones.size();
                // Copy ibm flat for back-compat consumers
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
            (void)vertCount;
        }
    }

    // Build images — tinygltf has already decoded each model.image via its
    // bundled stb_image into 8-bit RGBA when component == 4, or R/RG/RGB
    // variants otherwise. Normalize everything to RGBA8 for uniform consumers.
    scene.images.reserve(model.images.size());
    for (const auto& gi : model.images) {
        Image img;
        img.name = gi.name;
        img.mimeType = gi.mimeType;
        img.width  = gi.width;
        img.height = gi.height;

        if (gi.width > 0 && gi.height > 0 && !gi.image.empty()) {
            const size_t px = (size_t)gi.width * (size_t)gi.height;
            img.data.resize(px * 4);
            const int c = gi.component;
            const int bits = gi.bits;   // 8 or 16
            const uint8_t* src = gi.image.data();

            auto get8 = [&](size_t i, int channel) -> uint8_t {
                if (bits == 16) {
                    uint16_t v = ((const uint16_t*)src)[i * c + channel];
                    return (uint8_t)(v >> 8);
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

    // Build materials. Resolve each material's baseColorTexture to the
    // underlying image index (via model.textures[tex].source) so consumers
    // don't have to know about the glTF texture indirection.
    scene.materials.reserve(model.materials.size());
    for (const auto& gm : model.materials) {
        Material mat;
        mat.name = gm.name;
        const auto& pbr = gm.pbrMetallicRoughness;
        if (pbr.baseColorFactor.size() == 4) {
            mat.baseColorFactor[0] = (float)pbr.baseColorFactor[0];
            mat.baseColorFactor[1] = (float)pbr.baseColorFactor[1];
            mat.baseColorFactor[2] = (float)pbr.baseColorFactor[2];
            mat.baseColorFactor[3] = (float)pbr.baseColorFactor[3];
        }
        int texIdx = pbr.baseColorTexture.index;
        if (texIdx >= 0 && texIdx < (int)model.textures.size()) {
            int src = model.textures[texIdx].source;
            if (src >= 0 && src < (int)scene.images.size()) {
                mat.baseColorTexture = src;
            }
        }
        scene.materials.push_back(std::move(mat));
    }

    // Build animations.
    for (const auto& ganim : model.animations) {
        Animation anim;
        anim.name = ganim.name;
        anim.duration = 0.0f;
        int associatedSkel = -1;

        for (const auto& gch : ganim.channels) {
            AnimChannel ch;

            // Find which skeleton owns the target node.
            int boneIdx = -1;
            int skelIdx = -1;
            for (size_t si = 0; si < nodeToBone.size(); ++si) {
                auto it = nodeToBone[si].find(gch.target_node);
                if (it != nodeToBone[si].end()) {
                    boneIdx = it->second;
                    skelIdx = (int)si;
                    break;
                }
            }
            if (boneIdx < 0) continue; // not a skeletal channel
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
            ch.times.assign(tsrc, tsrc + accIn.count);
            if (!ch.times.empty()) {
                anim.duration = std::max(anim.duration, ch.times.back());
            }

            const auto& accOut = model.accessors[gs.output];
            const auto* vsrc = reinterpret_cast<const float*>(accessorData(model, gs.output));
            int stride = (ch.path == AnimChannel::Path::Rotation) ? 4 : 3;
            // CubicSpline packs (in-tangent, value, out-tangent) — 3x stride per key.
            int packing = (ch.interp == AnimChannel::Interp::CubicSpline) ? 3 : 1;
            ch.values.assign(vsrc, vsrc + accOut.count * stride * packing);

            anim.channels.push_back(std::move(ch));
        }

        scene.animations.push_back(std::move(anim));
        scene.animationSkeleton.push_back(associatedSkel);
    }

    return scene;
}

// ---- save helpers ----------------------------------------------------------

struct SaveCtx {
    tinygltf::Model model;
    tinygltf::Buffer buffer;

    int appendView(const void* data, size_t bytes, int target = 0) {
        size_t offset = buffer.data.size();
        buffer.data.resize(offset + bytes);
        std::memcpy(buffer.data.data() + offset, data, bytes);
        tinygltf::BufferView bv;
        bv.buffer = 0;
        bv.byteOffset = offset;
        bv.byteLength = bytes;
        if (target) bv.target = target;
        model.bufferViews.push_back(bv);
        return (int)model.bufferViews.size() - 1;
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
        return (int)model.accessors.size() - 1;
    }
};

static bool saveGLTFImpl(const MeshData& mesh,
                         const SkinData* skin,
                         const Skeleton* skeleton,
                         const std::vector<Animation>& animations,
                         const std::string& path) {
    if (mesh.empty()) return false;

    SaveCtx ctx;
    ctx.model.asset.version = "2.0";
    ctx.model.asset.generator = "bromesh";

    tinygltf::Mesh gltfMesh;
    tinygltf::Primitive prim;

    // --- positions (with min/max, required by spec) ---
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

    if (mesh.hasNormals()) {
        int v = ctx.appendView(mesh.normals.data(),
                               mesh.normals.size() * sizeof(float),
                               TINYGLTF_TARGET_ARRAY_BUFFER);
        int a = ctx.addAccessor(v, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_VEC3,
                                mesh.vertexCount());
        prim.attributes["NORMAL"] = a;
    }

    if (mesh.hasUVs()) {
        int v = ctx.appendView(mesh.uvs.data(),
                               mesh.uvs.size() * sizeof(float),
                               TINYGLTF_TARGET_ARRAY_BUFFER);
        int a = ctx.addAccessor(v, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_VEC2,
                                mesh.vertexCount());
        prim.attributes["TEXCOORD_0"] = a;
    }

    if (mesh.hasColors()) {
        int v = ctx.appendView(mesh.colors.data(),
                               mesh.colors.size() * sizeof(float),
                               TINYGLTF_TARGET_ARRAY_BUFFER);
        int a = ctx.addAccessor(v, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_VEC4,
                                mesh.vertexCount());
        prim.attributes["COLOR_0"] = a;
    }

    // --- skin attributes (JOINTS_0 as u16, WEIGHTS_0 as f32) ---
    const bool skinned = skin && skeleton && !skeleton->bones.empty()
                         && skin->boneIndices.size() == mesh.vertexCount() * 4;
    if (skinned) {
        std::vector<uint16_t> joints(skin->boneIndices.size());
        for (size_t i = 0; i < joints.size(); ++i)
            joints[i] = (uint16_t)skin->boneIndices[i];
        int v = ctx.appendView(joints.data(),
                               joints.size() * sizeof(uint16_t),
                               TINYGLTF_TARGET_ARRAY_BUFFER);
        int a = ctx.addAccessor(v, TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT,
                                TINYGLTF_TYPE_VEC4, mesh.vertexCount());
        prim.attributes["JOINTS_0"] = a;

        int v2 = ctx.appendView(skin->boneWeights.data(),
                                skin->boneWeights.size() * sizeof(float),
                                TINYGLTF_TARGET_ARRAY_BUFFER);
        int a2 = ctx.addAccessor(v2, TINYGLTF_COMPONENT_TYPE_FLOAT,
                                 TINYGLTF_TYPE_VEC4, mesh.vertexCount());
        prim.attributes["WEIGHTS_0"] = a2;
    }

    // --- indices ---
    {
        int v = ctx.appendView(mesh.indices.data(),
                               mesh.indices.size() * sizeof(uint32_t),
                               TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER);
        int a = ctx.addAccessor(v, TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT,
                                TINYGLTF_TYPE_SCALAR, mesh.indices.size());
        prim.indices = a;
    }
    prim.mode = TINYGLTF_MODE_TRIANGLES;
    gltfMesh.primitives.push_back(prim);
    ctx.model.meshes.push_back(gltfMesh);

    // --- nodes ---
    // Mesh node is node 0 (no skin) or comes after skeleton nodes (with skin).
    int meshNodeIdx = 0;
    int firstBoneNode = -1;
    if (skinned) {
        // One node per bone
        firstBoneNode = (int)ctx.model.nodes.size();
        for (const auto& b : skeleton->bones) {
            tinygltf::Node n;
            n.name = b.name;
            n.translation = { b.localT[0], b.localT[1], b.localT[2] };
            n.rotation    = { b.localR[0], b.localR[1], b.localR[2], b.localR[3] };
            n.scale       = { b.localS[0], b.localS[1], b.localS[2] };
            ctx.model.nodes.push_back(n);
        }
        // Link children
        for (size_t i = 0; i < skeleton->bones.size(); ++i) {
            int parent = skeleton->bones[i].parent;
            if (parent >= 0) {
                ctx.model.nodes[firstBoneNode + parent].children.push_back(firstBoneNode + (int)i);
            }
        }

        // Inverse bind matrices accessor
        std::vector<float> ibm(skeleton->bones.size() * 16);
        for (size_t i = 0; i < skeleton->bones.size(); ++i)
            std::memcpy(&ibm[i * 16], skeleton->bones[i].inverseBind, 16 * sizeof(float));
        int vIbm = ctx.appendView(ibm.data(), ibm.size() * sizeof(float));
        int aIbm = ctx.addAccessor(vIbm, TINYGLTF_COMPONENT_TYPE_FLOAT,
                                    TINYGLTF_TYPE_MAT4, skeleton->bones.size());

        tinygltf::Skin gskin;
        gskin.inverseBindMatrices = aIbm;
        for (size_t i = 0; i < skeleton->bones.size(); ++i)
            gskin.joints.push_back(firstBoneNode + (int)i);
        ctx.model.skins.push_back(gskin);
    }

    // Mesh node
    {
        tinygltf::Node n;
        n.mesh = 0;
        if (skinned) n.skin = 0;
        ctx.model.nodes.push_back(n);
        meshNodeIdx = (int)ctx.model.nodes.size() - 1;
    }

    // Scene root nodes: mesh node + root bones
    tinygltf::Scene gltfScene;
    gltfScene.nodes.push_back(meshNodeIdx);
    if (skinned) {
        for (size_t i = 0; i < skeleton->bones.size(); ++i) {
            if (skeleton->bones[i].parent < 0) {
                gltfScene.nodes.push_back(firstBoneNode + (int)i);
            }
        }
    }
    ctx.model.scenes.push_back(gltfScene);
    ctx.model.defaultScene = 0;

    // --- animations ---
    if (skinned) {
        for (const auto& anim : animations) {
            tinygltf::Animation gAnim;
            gAnim.name = anim.name;
            for (const auto& ch : anim.channels) {
                if (ch.boneIndex < 0 || ch.boneIndex >= (int)skeleton->bones.size()) continue;

                int vIn = ctx.appendView(ch.times.data(),
                                         ch.times.size() * sizeof(float));
                // For min/max on input accessor (required by spec):
                float tmin = ch.times.empty() ? 0.f : ch.times.front();
                float tmax = ch.times.empty() ? 0.f : ch.times.back();
                int aIn = ctx.addAccessor(vIn, TINYGLTF_COMPONENT_TYPE_FLOAT,
                                          TINYGLTF_TYPE_SCALAR, ch.times.size(),
                                          {(double)tmin}, {(double)tmax});

                int vOut = ctx.appendView(ch.values.data(),
                                          ch.values.size() * sizeof(float));
                int outType = (ch.path == AnimChannel::Path::Rotation)
                                  ? TINYGLTF_TYPE_VEC4 : TINYGLTF_TYPE_VEC3;
                int stride = (ch.path == AnimChannel::Path::Rotation) ? 4 : 3;
                int packing = (ch.interp == AnimChannel::Interp::CubicSpline) ? 3 : 1;
                size_t outCount = (stride == 0) ? 0 : (ch.values.size() / (stride * packing));
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
                gCh.sampler = (int)gAnim.samplers.size() - 1;
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
    }

    ctx.model.buffers.push_back(std::move(ctx.buffer));

    tinygltf::TinyGLTF writer;
    bool binary = path.size() >= 4 && path.substr(path.size() - 4) == ".glb";
    return writer.WriteGltfSceneToFile(&ctx.model, path, true, true, true, binary);
}

#endif // BROMESH_HAS_GLTF

// ---- public API ------------------------------------------------------------

#if !BROMESH_HAS_GLTF
GltfScene loadGLTF(const std::string&) { return {}; }
bool saveGLTF(const MeshData&, const std::string&) { return false; }
bool saveGLTF(const MeshData&, const SkinData*, const Skeleton*,
              const std::vector<Animation>&, const std::string&) { return false; }
#else
bool saveGLTF(const MeshData& mesh, const std::string& path) {
    return saveGLTFImpl(mesh, nullptr, nullptr, {}, path);
}
bool saveGLTF(const MeshData& mesh,
              const SkinData* skin,
              const Skeleton* skeleton,
              const std::vector<Animation>& animations,
              const std::string& path) {
    return saveGLTFImpl(mesh, skin, skeleton, animations, path);
}
#endif

} // namespace bromesh
