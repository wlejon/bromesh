#include "host_mesh_internal.h"

namespace bromesh::api {

HostClass g_poseClass;
HostClass g_animationClass;

namespace {

Value animChannelsToObject(const bromesh::Animation& anim) {
    return hostArrayOf(anim.channels.size(), [&](size_t i) {
        const auto& ch = anim.channels[i];
        ObjectBuilder obj;
        obj.set("boneIndex", static_cast<double>(ch.boneIndex));
        obj.set("path", ch.path == bromesh::AnimChannel::Path::Translation ? "translation" :
                        ch.path == bromesh::AnimChannel::Path::Rotation ? "rotation" : "scale");
        obj.set("interp", ch.interp == bromesh::AnimChannel::Interp::Step ? "STEP" :
                          ch.interp == bromesh::AnimChannel::Interp::CubicSpline ? "CUBICSPLINE" : "LINEAR");
        ev::Persistent t(makeFloat32Array(ch.times.data(), ch.times.size()));
        obj.set("times", t.get());
        ev::Persistent v(makeFloat32Array(ch.values.data(), ch.values.size()));
        obj.set("values", v.get());
        return obj.build();
    });
}

} // namespace

void initRiggingAnim(HostClass& poseCls, HostClass& animCls, HostClass& meshCls,
                     ObjectBuilder& ikBuilder) {
    // =========================================================================
    // Pose Class
    // =========================================================================
    poseCls.install("Pose", 2, [](Value, std::span<const Value> a) -> Value {
        auto h = std::make_unique<HostPose>();
        if (!a.empty()) {
            if (auto* skel = unwrapSkeleton(a[0])) {
                h->pose = bromesh::bindPose(skel->skeleton);
            } else if (ev::isNumber(a[0])) {
                size_t bc = static_cast<size_t>(ev::toDouble(a[0]));
                h->pose.data.assign(bc * 10, 0.0f);
                for (size_t i = 0; i < bc; ++i) {
                    h->pose.data[i * 10 + 6] = 1.0f; // rw = 1
                    h->pose.data[i * 10 + 7] = 1.0f; // sx = 1
                    h->pose.data[i * 10 + 8] = 1.0f; // sy = 1
                    h->pose.data[i * 10 + 9] = 1.0f; // sz = 1
                }
            } else {
                std::vector<float> d = toFloatVector(a[0]);
                if (a.size() > 1 && ev::isNumber(a[1])) {
                    size_t bc = static_cast<size_t>(ev::toDouble(a[1]));
                    if (d.size() < bc * 10) d.resize(bc * 10, 0.0f);
                }
                h->pose.data = std::move(d);
            }
        }
        return g_poseClass.createInstance(std::move(h));
    }, [](ObjectBuilder& proto) {
        proto.accessor("data",
            [](Value self, std::span<const Value>) -> Value {
                auto* p = unwrapPose(self);
                if (!p) return ev::undefined();
                return makeFloat32Array(p->pose.data.data(), p->pose.data.size());
            },
            [](Value self, std::span<const Value> a) -> Value {
                auto* p = unwrapPose(self);
                if (!p || a.empty()) return ev::undefined();
                p->pose.data = toFloatVector(a[0]);
                return ev::undefined();
            });

        proto.accessor("boneCount", [](Value self, std::span<const Value>) -> Value {
            auto* p = unwrapPose(self);
            return ev::fromDouble(p ? static_cast<double>(p->pose.boneCount()) : 0.0);
        });

        proto.def("computeWorldMatrices", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* p = unwrapPose(self);
            if (!p) return ev::throwTypeError("Pose.computeWorldMatrices: not a Pose instance");
            if (a.empty()) return ev::throwTypeError("Pose.computeWorldMatrices: Skeleton required");
            auto* s = unwrapSkeleton(a[0]);
            if (!s) return ev::throwTypeError("Pose.computeWorldMatrices: argument must be a Skeleton");
            std::vector<float> outWorld;
            bromesh::computeWorldMatrices(s->skeleton, p->pose, outWorld);
            return makeFloat32Array(outWorld.data(), outWorld.size());
        });

        proto.def("computeSkinningMatrices", 1, [](Value self, std::span<const Value> a) -> Value {
            auto* p = unwrapPose(self);
            if (!p) return ev::throwTypeError("Pose.computeSkinningMatrices: not a Pose instance");
            if (a.empty()) return ev::throwTypeError("Pose.computeSkinningMatrices: Skeleton required");
            auto* s = unwrapSkeleton(a[0]);
            if (!s) return ev::throwTypeError("Pose.computeSkinningMatrices: argument must be a Skeleton");
            std::vector<float> outSkinning;
            bromesh::computeSkinningMatrices(s->skeleton, p->pose, outSkinning);
            return makeFloat32Array(outSkinning.data(), outSkinning.size());
        });

        proto.def("socketWorld", 2, [](Value self, std::span<const Value> a) -> Value {
            auto* p = unwrapPose(self);
            if (!p || a.size() < 2) return ev::null();
            auto* s = unwrapSkeleton(a[0]);
            if (!s) return ev::null();
            std::string name = ev::toUtf8(a[1]);
            auto res = bromesh::socketWorldMatrix(s->skeleton, p->pose, name);
            if (!res.has_value()) return ev::null();
            return makeFloat32Array(res->data(), 16);
        });

        proto.def("clone", 0, [](Value self, std::span<const Value>) -> Value {
            auto* p = unwrapPose(self);
            if (!p) return ev::throwTypeError("Pose.clone: not a Pose instance");
            return wrapPose(p->pose);
        });
    });

    // Static blend on Pose
    poseCls.setStatic("blend", ev::makeFunction([](Value, std::span<const Value> a) -> Value {
        if (a.size() < 3) return ev::throwTypeError("Pose.blend: a, b, weight required");
        auto* pa = unwrapPose(a[0]);
        auto* pb = unwrapPose(a[1]);
        if (!pa || !pb) return ev::throwTypeError("Pose.blend: arguments must be Pose instances");
        float w = static_cast<float>(numAt(a, 2));
        const uint8_t* mask = nullptr;
        std::vector<uint8_t> maskVec;
        if (a.size() > 3 && !ev::isUndefined(a[3])) {
            maskVec = toUint8Vector(a[3]);
            if (!maskVec.empty()) mask = maskVec.data();
        }
        bromesh::Pose res = pa->pose;
        bromesh::blendPoses(res, pb->pose, w, mask);
        return wrapPose(std::move(res));
    }, 4, "blend"));

    poseCls.setStatic("blendN", ev::makeFunction([](Value, std::span<const Value> a) -> Value {
        if (a.size() < 2) return ev::throwTypeError("Pose.blendN: poses and weights required");
        std::vector<const bromesh::Pose*> posePtrs;
        if (ev::isObject(a[0])) {
            Value lenVal = ev::getProperty(a[0], "length");
            size_t n = ev::isNumber(lenVal) ? static_cast<size_t>(ev::toDouble(lenVal)) : 0;
            for (size_t i = 0; i < n; ++i) {
                Value elem = ev::getElement(a[0], static_cast<uint32_t>(i));
                auto* p = unwrapPose(elem);
                if (p) posePtrs.push_back(&p->pose);
            }
        }
        if (posePtrs.empty()) return ev::throwTypeError("Pose.blendN: no valid poses provided");

        std::vector<float> weights = toFloatVector(a[1]);
        if (weights.size() != posePtrs.size()) {
            return ev::throwTypeError("Pose.blendN: weights count must match poses count");
        }

        const uint8_t* mask = nullptr;
        std::vector<uint8_t> maskVec;
        if (a.size() > 2 && !ev::isUndefined(a[2])) {
            maskVec = toUint8Vector(a[2]);
            if (!maskVec.empty()) mask = maskVec.data();
        }

        bromesh::Pose out;
        bromesh::blendPosesN(posePtrs.data(), weights.data(), posePtrs.size(), out, mask);
        return wrapPose(std::move(out));
    }, 3, "blendN"));

    // =========================================================================
    // AnimationClip Class
    // =========================================================================
    animCls.install("AnimationClip", 1, [](Value, std::span<const Value> a) -> Value {
        auto h = std::make_unique<HostAnimation>();
        if (!a.empty() && ev::isObject(a[0])) {
            Value opts = a[0];
            Value nameVal = ev::getProperty(opts, "name");
            if (ev::isString(nameVal)) h->animation.name = ev::toUtf8(nameVal);
            Value durVal = ev::getProperty(opts, "duration");
            if (ev::isNumber(durVal)) h->animation.duration = static_cast<float>(ev::toDouble(durVal));

            Value chsVal = ev::getProperty(opts, "channels");
            if (ev::isObject(chsVal)) {
                Value lenVal = ev::getProperty(chsVal, "length");
                size_t n = ev::isNumber(lenVal) ? static_cast<size_t>(ev::toDouble(lenVal)) : 0;
                for (size_t i = 0; i < n; ++i) {
                    Value chVal = ev::getElement(chsVal, static_cast<uint32_t>(i));
                    if (ev::isObject(chVal)) {
                        bromesh::AnimChannel ch;
                        Value biVal = ev::getProperty(chVal, "boneIndex");
                        if (ev::isNumber(biVal)) ch.boneIndex = static_cast<int>(ev::toDouble(biVal));
                        Value pVal = ev::getProperty(chVal, "path");
                        if (ev::isString(pVal)) {
                            std::string ps = ev::toUtf8(pVal);
                            if (ps == "rotation") ch.path = bromesh::AnimChannel::Path::Rotation;
                            else if (ps == "scale") ch.path = bromesh::AnimChannel::Path::Scale;
                            else ch.path = bromesh::AnimChannel::Path::Translation;
                        }
                        Value tVal = ev::getProperty(chVal, "times");
                        if (!ev::isUndefined(tVal)) ch.times = toFloatVector(tVal);
                        Value vVal = ev::getProperty(chVal, "values");
                        if (!ev::isUndefined(vVal)) ch.values = toFloatVector(vVal);
                        h->animation.channels.push_back(std::move(ch));
                    }
                }
            }
        }
        return g_animationClass.createInstance(std::move(h));
    }, [](ObjectBuilder& proto) {
        proto.accessor("name",
            [](Value self, std::span<const Value>) -> Value {
                auto* a = unwrapAnimation(self);
                return ev::fromUtf8(a ? a->animation.name : "");
            },
            [](Value self, std::span<const Value> v) -> Value {
                auto* a = unwrapAnimation(self);
                if (a && !v.empty()) a->animation.name = ev::toUtf8(v[0]);
                return ev::undefined();
            });

        proto.accessor("duration",
            [](Value self, std::span<const Value>) -> Value {
                auto* a = unwrapAnimation(self);
                return ev::fromDouble(a ? static_cast<double>(a->animation.duration) : 0.0);
            },
            [](Value self, std::span<const Value> v) -> Value {
                auto* a = unwrapAnimation(self);
                if (a && !v.empty()) a->animation.duration = static_cast<float>(numAt(v, 0));
                return ev::undefined();
            });

        proto.accessor("channelCount", [](Value self, std::span<const Value>) -> Value {
            auto* a = unwrapAnimation(self);
            return ev::fromDouble(a ? static_cast<double>(a->animation.channels.size()) : 0.0);
        });

        proto.accessor("channels", [](Value self, std::span<const Value>) -> Value {
            auto* a = unwrapAnimation(self);
            if (!a) return hostArrayOf(std::span<const Value>{});
            return animChannelsToObject(a->animation);
        });

        proto.def("evaluate", 3, [](Value self, std::span<const Value> a) -> Value {
            auto* anim = unwrapAnimation(self);
            if (!anim) return ev::throwTypeError("AnimationClip.evaluate: not an instance");
            if (a.empty()) return ev::throwTypeError("AnimationClip.evaluate: skeleton required");
            auto* skel = unwrapSkeleton(a[0]);
            if (!skel) return ev::throwTypeError("AnimationClip.evaluate: argument must be a Skeleton");
            float t = a.size() > 1 ? static_cast<float>(numAt(a, 1)) : 0.0f;
            bool loop = a.size() > 2 ? boolAt(a, 2) : true;
            bromesh::Pose p = bromesh::evaluateAnimation(skel->skeleton, anim->animation, t, loop);
            return wrapPose(std::move(p));
        });
    });

    animCls.alias("Animation");
    animCls.alias("SkeletalAnimation");

    // =========================================================================
    // IK Solvers (IK namespace)
    // =========================================================================
    ikBuilder.def("twoBone", 7, [](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("IK.twoBone requires options or arguments");
        if (ev::isObject(a[0])) {
            Value opts = a[0];
            Value rootVal = ev::getProperty(opts, "rootPos");
            Value midVal = ev::getProperty(opts, "midPos");
            Value endVal = ev::getProperty(opts, "endPos");
            Value tgtVal = ev::getProperty(opts, "targetPos");
            Value poleVal = ev::getProperty(opts, "poleVector");

            auto r = toFloatVector(rootVal);
            auto m = toFloatVector(midVal);
            auto e = toFloatVector(endVal);
            auto t = toFloatVector(tgtVal);

            if (r.size() >= 3 && m.size() >= 3 && e.size() >= 3 && t.size() >= 3) {
                // Two-bone analytic triangle calculation
                bromath::Vec3 root{r[0], r[1], r[2]};
                bromath::Vec3 mid{m[0], m[1], m[2]};
                bromath::Vec3 end{e[0], e[1], e[2]};
                bromath::Vec3 target{t[0], t[1], t[2]};

                float l1 = bromath::vlen(mid - root);
                float l2 = bromath::vlen(end - mid);
                float d = bromath::vlen(target - root);
                if (d > l1 + l2 - 1e-4f) d = l1 + l2 - 1e-4f;
                if (d < std::fabs(l1 - l2) + 1e-4f) d = std::fabs(l1 - l2) + 1e-4f;

                float cosAngle = (l1 * l1 + d * d - l2 * l2) / (2.0f * l1 * d);
                cosAngle = std::clamp(cosAngle, -1.0f, 1.0f);
                float alpha = std::acos(cosAngle);

                bromath::Vec3 dir = bromath::vnorm(target - root);
                bromath::Vec3 pole{0, 1, 0};
                if (!ev::isUndefined(poleVal)) {
                    auto pv = toFloatVector(poleVal);
                    if (pv.size() >= 3) pole = {pv[0], pv[1], pv[2]};
                }
                bromath::Vec3 bendDir = bromath::vnorm(pole - bromath::vdot(pole, dir) * dir);
                bromath::Vec3 newMid = root + dir * (std::cos(alpha) * l1) + bendDir * (std::sin(alpha) * l1);

                ObjectBuilder res;
                {
                    const float midArr[3] = {newMid.x, newMid.y, newMid.z};
                    ev::Persistent p(makeFloat32Array(midArr, 3));
                    res.set("midPos", p.get());
                    res.set("endPos", tgtVal);
                }
                return res.build();
            }
        }
        return ev::undefined();
    });

    ikBuilder.def("solveTwoBone", 1, [](Value self, std::span<const Value> a) -> Value {
        return ev::call(ev::getProperty(self, "twoBone"), self, a).value;
    });

    ikBuilder.def("solveFabrik", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty() || !ev::isObject(a[0])) return hostArrayOf(std::span<const Value>{});
        Value opts = a[0];
        Value ptsVal = ev::getProperty(opts, "jointPositions");
        Value tgtVal = ev::getProperty(opts, "targetPos");
        auto tgt = toFloatVector(tgtVal);
        if (tgt.size() < 3 || !ev::isObject(ptsVal)) return ptsVal;

        Value lenVal = ev::getProperty(ptsVal, "length");
        size_t n = ev::isNumber(lenVal) ? static_cast<size_t>(ev::toDouble(lenVal)) : 0;
        std::vector<bromath::Vec3> points;
        points.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            auto pt = toFloatVector(ev::getElement(ptsVal, static_cast<uint32_t>(i)));
            if (pt.size() >= 3) points.push_back({pt[0], pt[1], pt[2]});
        }
        if (points.size() < 2) return ptsVal;

        bromath::Vec3 target{tgt[0], tgt[1], tgt[2]};
        std::vector<float> dists(points.size() - 1);
        for (size_t i = 0; i < points.size() - 1; ++i) {
            dists[i] = bromath::vlen(points[i + 1] - points[i]);
        }

        bromath::Vec3 origin = points[0];
        for (int iter = 0; iter < 10; ++iter) {
            // Backward
            points.back() = target;
            for (int i = static_cast<int>(points.size()) - 2; i >= 0; --i) {
                bromath::Vec3 r = bromath::vnorm(points[i] - points[i + 1]);
                points[i] = points[i + 1] + r * dists[i];
            }
            // Forward
            points[0] = origin;
            for (size_t i = 0; i < points.size() - 1; ++i) {
                bromath::Vec3 r = bromath::vnorm(points[i + 1] - points[i]);
                points[i + 1] = points[i] + r * dists[i];
            }
            if (bromath::vlen(points.back() - target) < 1e-3f) break;
        }

        return hostArrayOf(points.size(), [&](size_t i) {
            const float pt[3] = {points[i].x, points[i].y, points[i].z};
            return makeFloat32Array(pt, 3);
        });
    });

    ikBuilder.def("fabrik", 1, [](Value self, std::span<const Value> a) -> Value {
        return ev::call(ev::getProperty(self, "solveFabrik"), self, a).value;
    });

    ikBuilder.def("solveLookAt", 1, [](Value, std::span<const Value> a) -> Value {
        if (a.empty() || !ev::isObject(a[0])) return hostArrayOf(std::span<const Value>{});
        Value opts = a[0];
        auto head = toFloatVector(ev::getProperty(opts, "headPos"));
        auto tgt = toFloatVector(ev::getProperty(opts, "targetPos"));
        if (head.size() < 3 || tgt.size() < 3) return hostArrayOf(std::span<const Value>{});

        bromath::Vec3 forward = bromath::vnorm(bromath::Vec3{tgt[0] - head[0], tgt[1] - head[1], tgt[2] - head[2]});
        bromath::Vec3 up{0, 1, 0};
        bromath::Vec3 right = bromath::vnorm(bromath::vcross(up, forward));
        up = bromath::vcross(forward, right);

        float rot[4] = {0, 0, 0, 1}; // Quaternion approximation
        float tr = right.x + up.y + forward.z;
        if (tr > 0) {
            float s = 0.5f / std::sqrt(tr + 1.0f);
            rot[3] = 0.25f / s;
            rot[0] = (up.z - forward.y) * s;
            rot[1] = (forward.x - right.z) * s;
            rot[2] = (right.y - up.x) * s;
        }
        return makeFloat32Array(rot, 4);
    });

    ikBuilder.def("lookAt", 1, [](Value self, std::span<const Value> a) -> Value {
        return ev::call(ev::getProperty(self, "solveLookAt"), self, a).value;
    });

    // =========================================================================
    // Retargeting
    // =========================================================================
    animCls.setStatic("retarget", ev::makeFunction([](Value, std::span<const Value> a) -> Value {
        if (a.size() < 3) return ev::throwTypeError("AnimationClip.retarget: anim, srcSkel, dstSkel required");
        auto* anim = unwrapAnimation(a[0]);
        auto* src = unwrapSkeleton(a[1]);
        auto* dst = unwrapSkeleton(a[2]);
        if (!anim || !src || !dst) return ev::throwTypeError("AnimationClip.retarget: invalid arguments");
        auto retargeted = bromesh::retargetAnimation(anim->animation, src->skeleton, dst->skeleton);
        return wrapAnimation(std::move(retargeted));
    }, 3, "retarget"));

    // =========================================================================
    // Mesh Rigging Extensions
    // =========================================================================
    meshCls.prototype(); // ensures proto exists
    ObjectBuilder meshProto(g_meshClass.prototype());

    meshProto.def("applySkinning", 2, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.applySkinning: not a Mesh instance");
        if (a.size() < 2) return ev::throwTypeError("Mesh.applySkinning: skinData and poseMatrices required");
        auto* s = unwrapSkinData(a[0]);
        if (!s) return ev::throwTypeError("Mesh.applySkinning: first argument must be SkinData");
        std::vector<float> mats = toFloatVector(a[1]);
        if (mats.size() < s->skin.boneCount * 16) {
            return ev::throwTypeError("Mesh.applySkinning: not enough matrix floats for bones");
        }
        bromesh::applySkinning(m->mesh, s->skin, mats.data());
        return self;
    });

    meshProto.def("applyMorphTarget", 4, [](Value self, std::span<const Value> a) -> Value {
        auto* m = unwrapMesh(self);
        if (!m) return ev::throwTypeError("Mesh.applyMorphTarget: not a Mesh instance");
        if (a.size() < 3) return ev::throwTypeError("Mesh.applyMorphTarget: name, deltaPositions, weight required");
        bromesh::MorphTarget mt;
        mt.name = ev::toUtf8(a[0]);
        mt.deltaPositions = toFloatVector(a[1]);
        if (a.size() > 3) {
            mt.deltaNormals = toFloatVector(a[2]);
            float w = static_cast<float>(numAt(a, 3));
            bromesh::applyMorphTarget(m->mesh, mt, w);
        } else {
            float w = static_cast<float>(numAt(a, 2));
            bromesh::applyMorphTarget(m->mesh, mt, w);
        }
        return self;
    });

#if BROMESH_HAS_GLTF
    meshCls.setStatic("loadGLTF", ev::makeFunction([](Value, std::span<const Value> a) -> Value {
        if (a.empty()) return ev::throwTypeError("Mesh.loadGLTF: path required");
        std::string path = ev::toUtf8(a[0]);
        auto scene = bromesh::loadGLTF(path);
        ObjectBuilder res;
        Value meshesArr = ev::createArray();
        for (size_t i = 0; i < scene.meshes.size(); ++i) {
            ev::Persistent item(wrapMesh(std::move(scene.meshes[i])));
            ev::setElement(meshesArr, static_cast<uint32_t>(i), item.get());
        }
        res.set("meshes", meshesArr);

        Value skinsArr = ev::createArray();
        for (size_t i = 0; i < scene.skins.size(); ++i) {
            ev::Persistent item(wrapSkinData(std::move(scene.skins[i])));
            ev::setElement(skinsArr, static_cast<uint32_t>(i), item.get());
        }
        res.set("skins", skinsArr);

        Value skelsArr = ev::createArray();
        for (size_t i = 0; i < scene.skeletons.size(); ++i) {
            ev::Persistent item(wrapSkeleton(std::move(scene.skeletons[i])));
            ev::setElement(skelsArr, static_cast<uint32_t>(i), item.get());
        }
        res.set("skeletons", skelsArr);

        Value animsArr = ev::createArray();
        for (size_t i = 0; i < scene.animations.size(); ++i) {
            ev::Persistent item(wrapAnimation(std::move(scene.animations[i])));
            ev::setElement(animsArr, static_cast<uint32_t>(i), item.get());
        }
        res.set("animations", animsArr);

        return res.build();
    }, 1, "loadGLTF"));
#endif
}

void ensureRiggingClassesInstalled() {
    static bool installed = false;
    if (installed) return;
    installed = true;

    ensureMeshClassesInstalled();

    ObjectBuilder ikObj;
    initRiggingCore(g_skinDataClass, g_skeletonClass, g_jointClass, g_skeletonRigClass, g_voxelChunkClass);
    initRiggingAnim(g_poseClass, g_animationClass, g_meshClass, ikObj);

    // Register IK on globalThis
    ev::registerGlobal("IK", ikObj.get());
    ev::GlobalValue gt = ev::globalValue("globalThis");
    if (gt.found && !gt.value.isUndefined() && ev::isObject(gt.value)) {
        ev::setProperty(gt.value, "IK", ikObj.get());
    }
}

} // namespace bromesh::api
