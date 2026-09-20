#include "host_mesh_internal.h"
#include <fstream>
#include <sstream>

namespace bromesh::api {

HostClass g_skinDataClass;
HostClass g_skeletonClass;
HostClass g_jointClass;
HostClass g_skeletonRigClass;
HostClass g_voxelChunkClass;

namespace {

Value makeSkinValidationObject(const bromesh::SkinValidation& sv) {
    ObjectBuilder obj;
    obj.set("clean", sv.clean());
    obj.set("vertexCount", static_cast<double>(sv.vertexCount));
    obj.set("orphanCount", static_cast<double>(sv.orphanCount));
    obj.set("badSumCount", static_cast<double>(sv.badSumCount));
    obj.set("nanCount", static_cast<double>(sv.nanCount));
    obj.set("maxSumDeviation", static_cast<double>(sv.maxSumDeviation));
    obj.set("maxInfluencesObserved", static_cast<double>(sv.maxInfluencesObserved));
    return obj.build();
}

bromesh::Landmarks landmarksFromObject(Value v) {
    bromesh::Landmarks lm;
    if (!ev::isObject(v)) return lm;
    Value pts = ev::getProperty(v, "points");
    Value root = ev::isObject(pts) ? pts : v;

    Value keysVal = ev::getProperty(ev::globalValue("Object").value, "keys");
    if (ev::isFunction(keysVal)) {
        const Value args[1] = {root};
        auto res = ev::call(keysVal, ev::undefined(), std::span<const Value>(args, 1));
        if (!res.thrown && ev::isObject(res.value)) {
            Value lenVal = ev::getProperty(res.value, "length");
            size_t n = ev::isNumber(lenVal) ? static_cast<size_t>(ev::toDouble(lenVal)) : 0;
            for (size_t i = 0; i < n; ++i) {
                Value key = ev::getElement(res.value, static_cast<uint32_t>(i));
                std::string k = ev::toUtf8(key);
                Value ptVal = ev::getProperty(root, k);
                std::vector<float> pt = toFloatVector(ptVal);
                if (pt.size() >= 3) {
                    lm.set(k, pt[0], pt[1], pt[2]);
                }
            }
        }
    }
    return lm;
}

Value landmarksToObject(const bromesh::Landmarks& lm) {
    ObjectBuilder obj;
    ObjectBuilder pts;
    for (const auto& [name, p] : lm.points) {
        const float pt[3] = {p[0], p[1], p[2]};
        ev::Persistent v(makeFloat32Array(pt, 3));
        pts.set(name, v.get());
    }
    obj.set("points", pts.build());
    return obj.build();
}

} // namespace

void initRiggingCore(HostClass& skinCls, HostClass& skelCls, HostClass& jointCls,
                     HostClass& rigCls, HostClass& voxelCls) {
    // =========================================================================
    // SkinData Class
    // =========================================================================
    skinCls.install("SkinData", 1, [](Value, std::span<const Value> a) -> Value {
        auto h = std::make_unique<HostSkinData>();
        if (!a.empty() && ev::isObject(a[0])) {
            Value opts = a[0];
            Value bwVal = ev::getProperty(opts, "boneWeights");
            if (ev::isUndefined(bwVal)) bwVal = ev::getProperty(opts, "weights");
            Value biVal = ev::getProperty(opts, "boneIndices");
            if (ev::isUndefined(biVal)) biVal = ev::getProperty(opts, "indices");
            Value ibmVal = ev::getProperty(opts, "inverseBindMatrices");
            Value bcVal = ev::getProperty(opts, "boneCount");

            if (!ev::isUndefined(bwVal)) h->skin.boneWeights = toFloatVector(bwVal);
            if (!ev::isUndefined(biVal)) {
                std::vector<uint32_t> idx32 = toUint32Vector(biVal);
                h->skin.boneIndices.assign(idx32.begin(), idx32.end());
            }
            if (!ev::isUndefined(ibmVal)) h->skin.inverseBindMatrices = toFloatVector(ibmVal);
            if (ev::isNumber(bcVal)) {
                h->skin.boneCount = static_cast<size_t>(ev::toDouble(bcVal));
            } else if (!h->skin.inverseBindMatrices.empty()) {
                h->skin.boneCount = h->skin.inverseBindMatrices.size() / 16;
            }
        }
        return g_skinDataClass.createInstance(std::move(h));
    }, [](ObjectBuilder& proto) {
        proto.accessor("boneWeights", [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkinData(self);
            if (!s) return ev::undefined();
            return makeFloat32Array(s->skin.boneWeights.data(), s->skin.boneWeights.size());
        });
        proto.accessor("weights", [](Value self, std::span<const Value> a) -> Value {
            return ev::call(ev::getProperty(self, "boneWeights"), self, a).value;
        });

        proto.accessor("boneIndices", [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkinData(self);
            if (!s) return ev::undefined();
            std::vector<uint32_t> idx32(s->skin.boneIndices.begin(), s->skin.boneIndices.end());
            return makeUint32Array(idx32.data(), idx32.size());
        });
        proto.accessor("indices", [](Value self, std::span<const Value> a) -> Value {
            return ev::call(ev::getProperty(self, "boneIndices"), self, a).value;
        });

        proto.accessor("inverseBindMatrices", [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkinData(self);
            if (!s) return ev::undefined();
            return makeFloat32Array(s->skin.inverseBindMatrices.data(), s->skin.inverseBindMatrices.size());
        });

        proto.accessor("boneCount", [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkinData(self);
            return ev::fromDouble(s ? static_cast<double>(s->skin.boneCount) : 0.0);
        });

        proto.accessor("vertexCount", [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkinData(self);
            return ev::fromDouble(s ? static_cast<double>(s->skin.boneWeights.size() / 4) : 0.0);
        });

        proto.accessor("maxWeights", [](Value, std::span<const Value>) -> Value {
            return ev::fromDouble(4.0);
        });

        proto.def("normalize", 0, [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkinData(self);
            if (!s) return ev::throwTypeError("SkinData.normalize: not an instance");
            for (size_t v = 0; v < s->skin.boneWeights.size(); v += 4) {
                float sum = s->skin.boneWeights[v] + s->skin.boneWeights[v+1] +
                            s->skin.boneWeights[v+2] + s->skin.boneWeights[v+3];
                if (sum > 1e-6f) {
                    float inv = 1.0f / sum;
                    s->skin.boneWeights[v] *= inv;
                    s->skin.boneWeights[v+1] *= inv;
                    s->skin.boneWeights[v+2] *= inv;
                    s->skin.boneWeights[v+3] *= inv;
                }
            }
            return self;
        });

        proto.def("clone", 0, [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkinData(self);
            if (!s) return ev::throwTypeError("SkinData.clone: not an instance");
            return wrapSkinData(s->skin);
        });

        proto.def("validate", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* s = unwrapSkinData(self);
            if (!s) return ev::throwTypeError("SkinData.validate: not an instance");
            if (!a.empty()) {
                auto* m = unwrapMesh(a[0]);
                if (m && !m->mesh.empty()) {
                    auto v = bromesh::validateSkin(m->mesh, s->skin);
                    return makeSkinValidationObject(v);
                }
            }
            auto v = bromesh::validateSkin(s->skin);
            return makeSkinValidationObject(v);
        });
    });

    skinCls.setStatic("validate", ev::makeFunction([](Value, std::span<const Value> a) -> Value {
        bromesh::MeshData mesh;
        bromesh::SkinData skin;
        if (a.size() > 0) {
            auto* m = unwrapMesh(a[0]);
            if (m) mesh = m->mesh;
        }
        if (a.size() > 1) {
            auto* s = unwrapSkinData(a[1]);
            if (s) skin = s->skin;
        }
        auto v = bromesh::validateSkin(mesh, skin);
        return makeSkinValidationObject(v);
    }, 2, "validate"));

    skinCls.setStatic("transfer", ev::makeFunction([](Value, std::span<const Value> a) -> Value {
        if (a.size() < 3) return ev::throwTypeError("SkinData.transfer: targetMesh, sourceMesh, sourceSkin required");
        auto* tgt = unwrapMesh(a[0]);
        auto* src = unwrapMesh(a[1]);
        auto* skin = unwrapSkinData(a[2]);
        if (!tgt || !src || !skin) return ev::throwTypeError("SkinData.transfer: invalid arguments");
        float maxDist = a.size() > 3 ? static_cast<float>(numAt(a, 3)) : 0.0f;
        auto res = bromesh::transferSkinWeights(tgt->mesh, src->mesh, skin->skin, maxDist);
        return wrapSkinData(res);
    }, 4, "transfer"));

    // =========================================================================
    // Joint Class
    // =========================================================================
    jointCls.install("Joint", 1, [](Value, std::span<const Value> a) -> Value {
        auto h = std::make_unique<HostJoint>();
        if (!a.empty() && ev::isObject(a[0])) {
            Value opts = a[0];
            Value nameVal = ev::getProperty(opts, "name");
            if (ev::isString(nameVal)) h->bone.name = ev::toUtf8(nameVal);
            Value parentVal = ev::getProperty(opts, "parent");
            if (ev::isNumber(parentVal)) h->bone.parent = static_cast<int>(ev::toDouble(parentVal));
            Value idxVal = ev::getProperty(opts, "index");
            if (ev::isNumber(idxVal)) h->index = static_cast<int>(ev::toDouble(idxVal));

            Value tVal = ev::getProperty(opts, "localT");
            if (ev::isUndefined(tVal)) tVal = ev::getProperty(opts, "translation");
            if (!ev::isUndefined(tVal)) {
                auto t = toFloatVector(tVal);
                if (t.size() >= 3) { h->bone.localT[0] = t[0]; h->bone.localT[1] = t[1]; h->bone.localT[2] = t[2]; }
            }

            Value rVal = ev::getProperty(opts, "localR");
            if (ev::isUndefined(rVal)) rVal = ev::getProperty(opts, "rotation");
            if (!ev::isUndefined(rVal)) {
                auto r = toFloatVector(rVal);
                if (r.size() >= 4) { h->bone.localR[0] = r[0]; h->bone.localR[1] = r[1]; h->bone.localR[2] = r[2]; h->bone.localR[3] = r[3]; }
            }

            Value sVal = ev::getProperty(opts, "localS");
            if (ev::isUndefined(sVal)) sVal = ev::getProperty(opts, "scale");
            if (!ev::isUndefined(sVal)) {
                auto s = toFloatVector(sVal);
                if (s.size() >= 3) { h->bone.localS[0] = s[0]; h->bone.localS[1] = s[1]; h->bone.localS[2] = s[2]; }
            }

            Value ibmVal = ev::getProperty(opts, "inverseBind");
            if (ev::isUndefined(ibmVal)) ibmVal = ev::getProperty(opts, "inverseBindMatrix");
            if (!ev::isUndefined(ibmVal)) {
                auto ibm = toFloatVector(ibmVal);
                if (ibm.size() >= 16) {
                    for (int i = 0; i < 16; ++i) h->bone.inverseBind[i] = ibm[i];
                }
            }
        }
        return g_jointClass.createInstance(std::move(h));
    }, [](ObjectBuilder& proto) {
        proto.accessor("name",
            [](Value self, std::span<const Value>) -> Value {
                auto* j = unwrapJoint(self);
                return ev::fromUtf8(j ? j->bone.name : "");
            },
            [](Value self, std::span<const Value> a) -> Value {
                auto* j = unwrapJoint(self);
                if (j && !a.empty()) j->bone.name = ev::toUtf8(a[0]);
                return ev::undefined();
            });

        proto.accessor("parent",
            [](Value self, std::span<const Value>) -> Value {
                auto* j = unwrapJoint(self);
                return ev::fromDouble(j ? static_cast<double>(j->bone.parent) : -1.0);
            },
            [](Value self, std::span<const Value> a) -> Value {
                auto* j = unwrapJoint(self);
                if (j && !a.empty()) j->bone.parent = i32At(a, 0);
                return ev::undefined();
            });

        proto.accessor("index",
            [](Value self, std::span<const Value>) -> Value {
                auto* j = unwrapJoint(self);
                return ev::fromDouble(j ? static_cast<double>(j->index) : -1.0);
            },
            [](Value self, std::span<const Value> a) -> Value {
                auto* j = unwrapJoint(self);
                if (j && !a.empty()) j->index = i32At(a, 0);
                return ev::undefined();
            });

        proto.accessor("localT", [](Value self, std::span<const Value>) -> Value {
            auto* j = unwrapJoint(self);
            if (!j) return ev::undefined();
            return makeFloat32Array(j->bone.localT, 3);
        });

        proto.accessor("localR", [](Value self, std::span<const Value>) -> Value {
            auto* j = unwrapJoint(self);
            if (!j) return ev::undefined();
            return makeFloat32Array(j->bone.localR, 4);
        });

        proto.accessor("localS", [](Value self, std::span<const Value>) -> Value {
            auto* j = unwrapJoint(self);
            if (!j) return ev::undefined();
            return makeFloat32Array(j->bone.localS, 3);
        });

        proto.accessor("inverseBind", [](Value self, std::span<const Value>) -> Value {
            auto* j = unwrapJoint(self);
            if (!j) return ev::undefined();
            return makeFloat32Array(j->bone.inverseBind, 16);
        });
    });

    // =========================================================================
    // Skeleton Class
    // =========================================================================
    skelCls.install("Skeleton", 1, [](Value, std::span<const Value> a) -> Value {
        auto h = std::make_unique<HostSkeleton>();
        if (!a.empty() && ev::isObject(a[0])) {
            Value opts = a[0];
            Value bonesVal = ev::getProperty(opts, "bones");
            if (ev::isObject(bonesVal)) {
                Value lenVal = ev::getProperty(bonesVal, "length");
                size_t n = ev::isNumber(lenVal) ? static_cast<size_t>(ev::toDouble(lenVal)) : 0;
                h->skeleton.bones.reserve(n);
                for (size_t i = 0; i < n; ++i) {
                    Value bVal = ev::getElement(bonesVal, static_cast<uint32_t>(i));
                    bromesh::Bone b;
                    auto* j = unwrapJoint(bVal);
                    if (j) {
                        b = j->bone;
                    } else if (ev::isObject(bVal)) {
                        Value nameVal = ev::getProperty(bVal, "name");
                        if (ev::isString(nameVal)) b.name = ev::toUtf8(nameVal);
                        Value pVal = ev::getProperty(bVal, "parent");
                        if (ev::isNumber(pVal)) b.parent = static_cast<int>(ev::toDouble(pVal));

                        Value tVal = ev::getProperty(bVal, "localT");
                        if (ev::isUndefined(tVal)) tVal = ev::getProperty(bVal, "translation");
                        if (!ev::isUndefined(tVal)) {
                            auto t = toFloatVector(tVal);
                            if (t.size() >= 3) { b.localT[0] = t[0]; b.localT[1] = t[1]; b.localT[2] = t[2]; }
                        }
                        Value rVal = ev::getProperty(bVal, "localR");
                        if (ev::isUndefined(rVal)) rVal = ev::getProperty(bVal, "rotation");
                        if (!ev::isUndefined(rVal)) {
                            auto r = toFloatVector(rVal);
                            if (r.size() >= 4) { b.localR[0] = r[0]; b.localR[1] = r[1]; b.localR[2] = r[2]; b.localR[3] = r[3]; }
                        }
                        Value sVal = ev::getProperty(bVal, "localS");
                        if (ev::isUndefined(sVal)) sVal = ev::getProperty(bVal, "scale");
                        if (!ev::isUndefined(sVal)) {
                            auto s = toFloatVector(sVal);
                            if (s.size() >= 3) { b.localS[0] = s[0]; b.localS[1] = s[1]; b.localS[2] = s[2]; }
                        }
                        Value ibmVal = ev::getProperty(bVal, "inverseBind");
                        if (ev::isUndefined(ibmVal)) ibmVal = ev::getProperty(bVal, "inverseBindMatrix");
                        if (!ev::isUndefined(ibmVal)) {
                            auto ibm = toFloatVector(ibmVal);
                            if (ibm.size() >= 16) {
                                for (int k = 0; k < 16; ++k) b.inverseBind[k] = ibm[k];
                            }
                        }
                    }
                    h->skeleton.bones.push_back(b);
                }
            }

            Value socketsVal = ev::getProperty(opts, "sockets");
            if (ev::isObject(socketsVal)) {
                Value lenVal = ev::getProperty(socketsVal, "length");
                size_t n = ev::isNumber(lenVal) ? static_cast<size_t>(ev::toDouble(lenVal)) : 0;
                for (size_t i = 0; i < n; ++i) {
                    Value sVal = ev::getElement(socketsVal, static_cast<uint32_t>(i));
                    if (ev::isObject(sVal)) {
                        bromesh::Socket sock;
                        Value nameVal = ev::getProperty(sVal, "name");
                        if (ev::isString(nameVal)) sock.name = ev::toUtf8(nameVal);
                        Value boneVal = ev::getProperty(sVal, "boneIndex");
                        if (ev::isUndefined(boneVal)) boneVal = ev::getProperty(sVal, "bone");
                        if (ev::isNumber(boneVal)) sock.bone = static_cast<int>(ev::toDouble(boneVal));
                        Value offVal = ev::getProperty(sVal, "offset");
                        if (!ev::isUndefined(offVal)) {
                            auto off = toFloatVector(offVal);
                            if (off.size() >= 16) {
                                for (int k = 0; k < 16; ++k) sock.offset[k] = off[k];
                            }
                        }
                        h->skeleton.sockets.push_back(sock);
                    }
                }
            }
        }
        return g_skeletonClass.createInstance(std::move(h));
    }, [](ObjectBuilder& proto) {
        proto.accessor("boneCount", [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkeleton(self);
            return ev::fromDouble(s ? static_cast<double>(s->skeleton.bones.size()) : 0.0);
        });

        proto.accessor("socketCount", [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkeleton(self);
            return ev::fromDouble(s ? static_cast<double>(s->skeleton.sockets.size()) : 0.0);
        });

        proto.accessor("bones", [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkeleton(self);
            if (!s) return hostArrayOf(std::span<const Value>{});
            return hostArrayOf(s->skeleton.bones.size(), [&](size_t i) {
                return wrapJoint(s->skeleton.bones[i], static_cast<int>(i));
            });
        });

        proto.accessor("sockets", [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkeleton(self);
            if (!s) return hostArrayOf(std::span<const Value>{});
            return hostArrayOf(s->skeleton.sockets.size(), [&](size_t i) {
                const auto& sock = s->skeleton.sockets[i];
                ObjectBuilder obj;
                obj.set("name", sock.name);
                obj.set("bone", static_cast<double>(sock.bone));
                obj.set("boneIndex", static_cast<double>(sock.bone));
                ev::Persistent off(makeFloat32Array(sock.offset, 16));
                obj.set("offset", off.get());
                return obj.build();
            });
        });

        proto.def("findBone", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* s = unwrapSkeleton(self);
            if (!s || a.empty()) return ev::fromDouble(-1.0);
            std::string name = ev::toUtf8(a[0]);
            return ev::fromDouble(static_cast<double>(s->skeleton.findBone(name)));
        });

        proto.def("findBoneBySuffix", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* s = unwrapSkeleton(self);
            if (!s || a.empty()) return ev::fromDouble(-1.0);
            std::string suffix = ev::toUtf8(a[0]);
            for (size_t i = 0; i < s->skeleton.bones.size(); ++i) {
                const std::string& bname = s->skeleton.bones[i].name;
                if (bname.size() >= suffix.size() &&
                    bname.compare(bname.size() - suffix.size(), suffix.size(), suffix) == 0) {
                    return ev::fromDouble(static_cast<double>(i));
                }
            }
            return ev::fromDouble(-1.0);
        });

        proto.def("findSocket", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* s = unwrapSkeleton(self);
            if (!s || a.empty()) return ev::fromDouble(-1.0);
            std::string name = ev::toUtf8(a[0]);
            for (size_t i = 0; i < s->skeleton.sockets.size(); ++i) {
                if (s->skeleton.sockets[i].name == name) return ev::fromDouble(static_cast<double>(i));
            }
            return ev::fromDouble(-1.0);
        });

        proto.def("addSocket", 3, [](Value self, std::span<const Value> a) -> Value {
            auto* s = unwrapSkeleton(self);
            if (!s || a.empty()) return ev::fromDouble(-1.0);
            bromesh::Socket sock;
            if (ev::isObject(a[0]) && a.size() == 1) {
                Value nameVal = ev::getProperty(a[0], "name");
                Value bVal = ev::getProperty(a[0], "bone");
                if (ev::isUndefined(bVal)) bVal = ev::getProperty(a[0], "boneIndex");
                Value offVal = ev::getProperty(a[0], "offset");
                sock.name = ev::toUtf8(nameVal);
                sock.bone = ev::isNumber(bVal) ? static_cast<int>(ev::toDouble(bVal)) : 0;
                auto off = toFloatVector(offVal);
                if (off.size() >= 16) for (int i = 0; i < 16; ++i) sock.offset[i] = off[i];
            } else {
                sock.name = ev::toUtf8(a[0]);
                sock.bone = a.size() > 1 ? i32At(a, 1) : 0;
                if (a.size() > 2) {
                    auto off = toFloatVector(a[2]);
                    if (off.size() >= 16) for (int i = 0; i < 16; ++i) sock.offset[i] = off[i];
                }
            }
            s->skeleton.sockets.push_back(sock);
            return ev::fromDouble(static_cast<double>(s->skeleton.sockets.size() - 1));
        });

        // addRigifySockets() -> number — append the standard attachment
        // sockets (hands, feet, head, spine) for a Rigify / Mixamo-named
        // skeleton, matching "ORG-"/"DEF-"/bare/"mixamorig:" bone spellings.
        // Returns how many were added. Dropped by the bronze port
        // (bro docs/transition-drift.md H7).
        proto.def("addRigifySockets", 0, [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkeleton(self);
            if (!s) return ev::throwTypeError("Skeleton.addRigifySockets: not an instance");
            return ev::fromDouble(static_cast<double>(bromesh::addRigifySockets(s->skeleton)));
        });

        proto.def("bindPose", 0, [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkeleton(self);
            if (!s) return ev::throwTypeError("Skeleton.bindPose: not an instance");
            return wrapPose(bromesh::bindPose(s->skeleton));
        });

        proto.def("clone", 0, [](Value self, std::span<const Value>) -> Value {
            auto* s = unwrapSkeleton(self);
            if (!s) return ev::throwTypeError("Skeleton.clone: not an instance");
            return wrapSkeleton(s->skeleton);
        });
    });

    skelCls.setStatic("fromBones", ev::makeFunction([](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return wrapSkeleton(bromesh::Skeleton{});
        ObjectBuilder opts;
        opts.set("bones", a[0]);
        const Value args[1] = {opts.build()};
        return ev::call(g_skeletonClass.constructor(), ev::undefined(), std::span<const Value>(args, 1)).value;
    }, 1, "fromBones"));

    // =========================================================================
    // SkeletonRig (Rig / RigSpec) Class
    // =========================================================================
    rigCls.install("SkeletonRig", 1, [](Value, std::span<const Value> a) -> Value {
        auto h = std::make_unique<HostSkeletonRig>();
        if (!a.empty()) {
            if (ev::isString(a[0])) {
                std::string type = ev::toUtf8(a[0]);
                if (type == "humanoid") h->spec = bromesh::builtinHumanoidSpec();
                else if (type == "quadruped") h->spec = bromesh::builtinQuadrupedSpec();
                else h->spec.name = type;
            } else if (ev::isObject(a[0])) {
                Value tVal = ev::getProperty(a[0], "type");
                if (ev::isUndefined(tVal)) tVal = ev::getProperty(a[0], "name");
                if (ev::isString(tVal)) {
                    std::string t = ev::toUtf8(tVal);
                    if (t == "humanoid") h->spec = bromesh::builtinHumanoidSpec();
                    else if (t == "quadruped") h->spec = bromesh::builtinQuadrupedSpec();
                    else h->spec.name = t;
                }
            }
        } else {
            h->spec = bromesh::builtinHumanoidSpec();
        }
        return g_skeletonRigClass.createInstance(std::move(h));
    }, [](ObjectBuilder& proto) {
        proto.accessor("type", [](Value self, std::span<const Value>) -> Value {
            auto* r = unwrapSkeletonRig(self);
            return ev::fromUtf8(r ? r->spec.name : "");
        });
        proto.accessor("name", [](Value self, std::span<const Value>) -> Value {
            auto* r = unwrapSkeletonRig(self);
            return ev::fromUtf8(r ? r->spec.name : "");
        });
        proto.accessor("boneCount", [](Value self, std::span<const Value>) -> Value {
            auto* r = unwrapSkeletonRig(self);
            return ev::fromDouble(r ? static_cast<double>(r->spec.bones.size()) : 0.0);
        });
        proto.accessor("landmarkCount", [](Value self, std::span<const Value>) -> Value {
            auto* r = unwrapSkeletonRig(self);
            return ev::fromDouble(r ? static_cast<double>(r->spec.landmarks.size()) : 0.0);
        });
        proto.def("toJSON", 0, [](Value self, std::span<const Value>) -> Value {
            auto* r = unwrapSkeletonRig(self);
            return ev::fromUtf8(r ? bromesh::serializeRigSpecJSON(r->spec) : "");
        });
        proto.def("landmarkNames", 0, [](Value self, std::span<const Value>) -> Value {
            auto* r = unwrapSkeletonRig(self);
            if (!r) return hostArrayOf(std::span<const Value>{});
            return hostArrayOf(r->spec.landmarks.size(), [&](size_t i) {
                return ev::fromUtf8(r->spec.landmarks[i].name);
            });
        });
        proto.def("boneNames", 0, [](Value self, std::span<const Value>) -> Value {
            auto* r = unwrapSkeletonRig(self);
            if (!r) return hostArrayOf(std::span<const Value>{});
            return hostArrayOf(r->spec.bones.size(), [&](size_t i) {
                return ev::fromUtf8(r->spec.bones[i].name);
            });
        });
    });

    // Aliases on globalThis
    rigCls.alias("RigSpec");
    rigCls.alias("Rig");

    // Static Rig methods
    auto bindRigStatic = [&](const char* name, uint32_t arity, ev::NativeFn fn) {
        rigCls.setStatic(name, ev::makeFunction(std::move(fn), arity, name));
    };

    bindRigStatic("spec", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return wrapSkeletonRig(bromesh::RigSpec{});
        std::string name = ev::toUtf8(a[0]);
        return wrapSkeletonRig(bromesh::builtinRigSpec(name));
    });

    bindRigStatic("specFromJSON", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return wrapSkeletonRig(bromesh::RigSpec{});
        std::string json = ev::toUtf8(a[0]);
        return wrapSkeletonRig(bromesh::parseRigSpecJSON(json));
    });

    bindRigStatic("specFromFile", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Rig.specFromFile: path required");
        std::string path = ev::toUtf8(a[0]);
        return wrapSkeletonRig(bromesh::loadRigSpecFile(path));
    });

    bindRigStatic("detectHumanoid", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Rig.detectHumanoid: mesh required");
        auto* m = unwrapMesh(a[0]);
        if (!m) return ev::throwTypeError("Rig.detectHumanoid: mesh required");
        auto lm = bromesh::detectHumanoidLandmarks(m->mesh);
        return landmarksToObject(lm);
    });

    bindRigStatic("detectLandmarks", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Rig.detectLandmarks: mesh required");
        auto* m = unwrapMesh(a[0]);
        if (!m) return ev::throwTypeError("Rig.detectLandmarks: mesh required");
        auto lm = bromesh::detectHumanoidLandmarks(m->mesh);
        return landmarksToObject(lm);
    });

    bindRigStatic("detectQuadruped", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Rig.detectQuadruped: mesh required");
        auto* m = unwrapMesh(a[0]);
        if (!m) return ev::throwTypeError("Rig.detectQuadruped: mesh required");
        auto lm = bromesh::detectQuadrupedLandmarks(m->mesh);
        return landmarksToObject(lm);
    });

    bindRigStatic("missingLandmarks", 2, [](Value, std::span<const Value> a) -> Value {
        if (a.size() < 2) return hostArrayOf(std::span<const Value>{});
        auto* r = unwrapSkeletonRig(a[0]);
        bromesh::RigSpec spec = r ? r->spec : bromesh::builtinHumanoidSpec();
        bromesh::Landmarks lm = landmarksFromObject(a[1]);
        auto missing = bromesh::missingLandmarks(spec, lm);
        return hostArrayOf(missing.size(), [&](size_t i) {
            return ev::fromUtf8(missing[i]);
        });
    });

    bindRigStatic("fitSkeleton", 3, [](Value, std::span<const Value> a) -> Value {
        if (a.size() < 3) return ev::throwTypeError("Rig.fitSkeleton: spec, landmarks, and mesh required");
        auto* r = unwrapSkeletonRig(a[0]);
        bromesh::RigSpec spec = r ? r->spec : bromesh::builtinHumanoidSpec();
        bromesh::Landmarks lm = landmarksFromObject(a[1]);
        auto* m = unwrapMesh(a[2]);
        if (!m) return ev::throwTypeError("Rig.fitSkeleton: mesh required");
        return wrapSkeleton(bromesh::fitSkeleton(spec, lm, m->mesh));
    });

    bindRigStatic("autoRig", 4, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Rig.autoRig: mesh required");
        auto* m = unwrapMesh(a[0]);
        if (!m) return ev::throwTypeError("Rig.autoRig: mesh must be a Mesh instance");

        bromesh::RigSpec spec = bromesh::builtinHumanoidSpec();
        bromesh::Landmarks lm;
        bool hasLandmarks = false;
        Value opts = ev::undefined();

        if (a.size() >= 3 && unwrapSkeletonRig(a[1])) {
            spec = unwrapSkeletonRig(a[1])->spec;
            if (ev::isObject(a[2])) {
                lm = landmarksFromObject(a[2]);
                hasLandmarks = true;
            }
            if (a.size() > 3 && ev::isObject(a[3])) opts = a[3];
        } else if (a.size() > 1 && ev::isObject(a[1])) {
            opts = a[1];
            Value specVal = ev::getProperty(opts, "spec");
            if (auto* r = unwrapSkeletonRig(specVal)) spec = r->spec;
            else {
                Value typeVal = ev::getProperty(opts, "rigType");
                if (ev::isString(typeVal)) {
                    spec = bromesh::builtinRigSpec(ev::toUtf8(typeVal));
                    if (spec.name.empty()) spec.name = ev::toUtf8(typeVal);
                }
            }
            Value lmVal = ev::getProperty(opts, "landmarks");
            if (ev::isObject(lmVal)) {
                lm = landmarksFromObject(lmVal);
                hasLandmarks = true;
            }
        }

        if (!hasLandmarks) {
            if (spec.name == "quadruped") lm = bromesh::detectQuadrupedLandmarks(m->mesh);
            else lm = bromesh::detectHumanoidLandmarks(m->mesh);
        }

        bromesh::WeightingOptions wopts;
        if (ev::isObject(opts)) {
            Value mVal = ev::getProperty(opts, "method");
            if (ev::isString(mVal)) {
                wopts.method = bromesh::parseWeightingMethod(ev::toUtf8(mVal).c_str());
            }
            Value smVal = ev::getProperty(opts, "smoothIterations");
            if (ev::isNumber(smVal)) wopts.smoothIterations = static_cast<int>(ev::toDouble(smVal));
        }

        auto res = bromesh::autoRig(m->mesh, spec, lm, wopts);
        ObjectBuilder obj;
        {
            ev::Persistent skel(wrapSkeleton(std::move(res.skeleton)));
            obj.set("skeleton", skel.get());
            ev::Persistent skin(wrapSkinData(std::move(res.skin)));
            obj.set("skin", skin.get());
        }
        obj.set("methodUsed", res.methodUsed == bromesh::WeightingMethod::BBW ? "bbw" :
                              res.methodUsed == bromesh::WeightingMethod::BoneHeat ? "boneHeat" : "voxelBind");
        obj.set("missingLandmarks", hostArrayOf(res.missingLandmarks.size(), [&](size_t i) {
            return ev::fromUtf8(res.missingLandmarks[i]);
        }));
        obj.set("warnings", hostArrayOf(res.warnings.size(), [&](size_t i) {
            return ev::fromUtf8(res.warnings[i]);
        }));
        return obj.build();
    });

    bindRigStatic("generateLocomotionCycle", 3, [](Value, std::span<const Value> a) -> Value {
        if (a.size() < 2) return ev::throwTypeError("Rig.generateLocomotionCycle: skeleton and spec required");
        auto* s = unwrapSkeleton(a[0]);
        auto* r = unwrapSkeletonRig(a[1]);
        if (!s || !r) return ev::throwTypeError("Rig.generateLocomotionCycle: invalid arguments");
        bromesh::LocomotionParams params;
        if (a.size() > 2 && ev::isString(a[2])) {
            params.gait.name = ev::toUtf8(a[2]);
        }
        auto anim = bromesh::generateLocomotionCycle(s->skeleton, r->spec, params);
        return wrapAnimation(std::move(anim));
    });

    bindRigStatic("transferWeights", 3, [](Value, std::span<const Value> a) -> Value {
        return ev::call(ev::getProperty(g_skinDataClass.constructor(), "transfer"), ev::undefined(), a).value;
    });

    // =========================================================================
    // VoxelChunk Class
    // =========================================================================
    voxelCls.install("VoxelChunk", 4, [](Value, std::span<const Value> a) -> Value {
        ArgReader r(a);
        int dx = r.getInt(0, 16);
        int dy = r.getInt(1, 16);
        int dz = r.getInt(2, 16);
        float cs = static_cast<float>(r.getDouble(3, 1.0));
        auto h = std::make_unique<HostVoxelChunk>();
        h->chunk = std::make_unique<bromesh::VoxelChunk>(dx > 0 ? dx : 16, dy > 0 ? dy : 16, dz > 0 ? dz : 16, cs > 0.0f ? cs : 1.0f);
        return g_voxelChunkClass.createInstance(std::move(h));
    }, [](ObjectBuilder& proto) {
        proto.accessor("sizeX", [](Value self, std::span<const Value>) -> Value {
            auto* v = unwrapVoxelChunk(self);
            return ev::fromDouble(v && v->chunk ? static_cast<double>(v->chunk->sizeX()) : 0.0);
        });
        proto.accessor("sizeY", [](Value self, std::span<const Value>) -> Value {
            auto* v = unwrapVoxelChunk(self);
            return ev::fromDouble(v && v->chunk ? static_cast<double>(v->chunk->sizeY()) : 0.0);
        });
        proto.accessor("sizeZ", [](Value self, std::span<const Value>) -> Value {
            auto* v = unwrapVoxelChunk(self);
            return ev::fromDouble(v && v->chunk ? static_cast<double>(v->chunk->sizeZ()) : 0.0);
        });
        proto.accessor("cellSize", [](Value self, std::span<const Value>) -> Value {
            auto* v = unwrapVoxelChunk(self);
            return ev::fromDouble(v && v->chunk ? static_cast<double>(v->chunk->cellSize()) : 1.0);
        });
        proto.accessor("isDirty",
            [](Value self, std::span<const Value>) -> Value {
                auto* v = unwrapVoxelChunk(self);
                return ev::fromBool(v && v->chunk ? v->chunk->isDirty() : false);
            },
            [](Value self, std::span<const Value> a) -> Value {
                auto* v = unwrapVoxelChunk(self);
                if (v && v->chunk) {
                    if (!a.empty() && ev::toBool(a[0])) v->chunk->markDirty();
                    else v->chunk->clearDirty();
                }
                return ev::undefined();
            });

        proto.def("clearDirty", 0, [](Value self, std::span<const Value>) -> Value {
            auto* v = unwrapVoxelChunk(self);
            if (v && v->chunk) v->chunk->clearDirty();
            return self;
        });
        proto.def("markDirty", 0, [](Value self, std::span<const Value>) -> Value {
            auto* v = unwrapVoxelChunk(self);
            if (v && v->chunk) v->chunk->markDirty();
            return self;
        });

        proto.def("set", 4, [](Value self, std::span<const Value> a) -> Value {
            auto* v = unwrapVoxelChunk(self);
            if (!v || !v->chunk) return ev::throwTypeError("VoxelChunk.set: not an instance");
            ArgReader r(a);
            v->chunk->setVoxel(r.getInt(0, 0), r.getInt(1, 0), r.getInt(2, 0), static_cast<uint8_t>(r.getInt(3, 0)));
            return self;
        });
        proto.def("setVoxel", 4, [](Value self, std::span<const Value> a) -> Value {
            return ev::call(ev::getProperty(self, "set"), self, a).value;
        });

        proto.def("get", 3, [](Value self, std::span<const Value> a) -> Value {
            auto* v = unwrapVoxelChunk(self);
            if (!v || !v->chunk) return ev::throwTypeError("VoxelChunk.get: not an instance");
            ArgReader r(a);
            return ev::fromDouble(static_cast<double>(v->chunk->getVoxel(r.getInt(0, 0), r.getInt(1, 0), r.getInt(2, 0))));
        });
        proto.def("getVoxel", 3, [](Value self, std::span<const Value> a) -> Value {
            return ev::call(ev::getProperty(self, "get"), self, a).value;
        });

        proto.def("fill", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* v = unwrapVoxelChunk(self);
            if (!v || !v->chunk) return ev::throwTypeError("VoxelChunk.fill: not an instance");
            v->chunk->fill(static_cast<uint8_t>(i32At(a, 0)));
            return self;
        });

        proto.def("markDirty", 0, [](Value self, std::span<const Value>) -> Value {
            auto* v = unwrapVoxelChunk(self);
            if (v && v->chunk) v->chunk->markDirty();
            return self;
        });
        proto.def("clearDirty", 0, [](Value self, std::span<const Value>) -> Value {
            auto* v = unwrapVoxelChunk(self);
            if (v && v->chunk) v->chunk->clearDirty();
            return self;
        });

        proto.def("data", 0, [](Value self, std::span<const Value>) -> Value {
            auto* v = unwrapVoxelChunk(self);
            if (!v || !v->chunk) return ev::throwTypeError("VoxelChunk.data: not an instance");
            size_t sz = static_cast<size_t>(v->chunk->sizeX()) * v->chunk->sizeY() * v->chunk->sizeZ();
            return makeUint8Array(v->chunk->data(), sz);
        });

        proto.def("setData", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* v = unwrapVoxelChunk(self);
            if (!v || !v->chunk || a.empty()) return ev::throwTypeError("VoxelChunk.setData: not an instance or empty");
            std::vector<uint8_t> d = toUint8Vector(a[0]);
            size_t sz = static_cast<size_t>(v->chunk->sizeX()) * v->chunk->sizeY() * v->chunk->sizeZ();
            std::memcpy(v->chunk->data(), d.data(), std::min(sz, d.size()));
            v->chunk->markDirty();
            return self;
        });

        proto.def("buildMesh", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* v = unwrapVoxelChunk(self);
            if (!v || !v->chunk) return ev::throwTypeError("VoxelChunk.buildMesh: not an instance");
            const float* pal = nullptr;
            int count = 0;
            std::vector<float> palVec;
            if (!a.empty()) {
                palVec = toFloatVector(a[0]);
                if (!palVec.empty()) {
                    pal = palVec.data();
                    count = a.size() > 1 ? i32At(a, 1) : static_cast<int>(palVec.size() / 4);
                }
            }
            return wrapMesh(v->chunk->buildMesh(pal, count));
        });
        proto.def("toMesh", 0, [](Value self, std::span<const Value> a) -> Value {
            return ev::call(ev::getProperty(self, "buildMesh"), self, a).value;
        });
    });
}

} // namespace bromesh::api
