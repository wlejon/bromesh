#include "mir_lane_emitter.h"
#include "bromesh/isosurface/jit/mir_mesh_kernels.h"
#include "sdf_compiler_internal.h"

#include <cmath>
#include <vector>

namespace bromesh {

void evaluateSdfGridReference(const SdfGraph& graph, int dimX, int dimY, int dimZ,
                              const bromath::AABB3& bounds, std::vector<float>& outField,
                              int rootNode) {
    if (dimX <= 0 || dimY <= 0 || dimZ <= 0) {
        outField.clear();
        return;
    }
    const int root = (rootNode >= 0) ? rootNode : graph.root();
    const float ox = bounds.min.x;
    const float oy = bounds.min.y;
    const float oz = bounds.min.z;
    const float sx = (dimX > 1) ? (bounds.max.x - bounds.min.x) / static_cast<float>(dimX - 1) : 1.0f;
    const float sy = (dimY > 1) ? (bounds.max.y - bounds.min.y) / static_cast<float>(dimY - 1) : 1.0f;
    const float sz = (dimZ > 1) ? (bounds.max.z - bounds.min.z) / static_cast<float>(dimZ - 1) : 1.0f;

    const size_t total = static_cast<size_t>(dimX) * static_cast<size_t>(dimY) * static_cast<size_t>(dimZ);
    outField.resize(total);
    for (int z = 0; z < dimZ; ++z) {
        const float pz = std::fma(static_cast<float>(z), sz, oz);
        for (int y = 0; y < dimY; ++y) {
            const float py = std::fma(static_cast<float>(y), sy, oy);
            for (int x = 0; x < dimX; ++x) {
                const float px = std::fma(static_cast<float>(x), sx, ox);
                const size_t idx = (static_cast<size_t>(z) * static_cast<size_t>(dimY) +
                                    static_cast<size_t>(y)) * static_cast<size_t>(dimX) +
                                   static_cast<size_t>(x);
                outField[idx] = evalScalarNode(graph, root, px, py, pz);
            }
        }
    }
}

#if BROMESH_HAS_BRASS_JIT
namespace mir_detail {

using namespace brass;

Type MirLaneEmitter::f32Type() const noexcept {
    return (lanes_ == 8) ? Type::f32x8() : (lanes_ == 4) ? Type::f32x4() : Type::f32();
}

Type MirLaneEmitter::i32Type() const noexcept {
    return (lanes_ == 8) ? Type::i32x8() : (lanes_ == 4) ? Type::i32x4() : Type::i32();
}

Value* MirLaneEmitter::broadcastF32(Value* scalar) {
    return (lanes_ == 1) ? scalar : b_.build_vbroadcast(f32Type(), scalar);
}

Value* MirLaneEmitter::cF32(float v) {
    return broadcastF32(kb_.const_f32(v));
}

Value* MirLaneEmitter::cI32(int32_t v) {
    Value* s = kb_.const_i32(v);
    return (lanes_ == 1) ? s : b_.build_vbroadcast(i32Type(), s);
}

Value* MirLaneEmitter::add(Value* a, Value* b) { return (lanes_ == 1) ? b_.build_add(a, b) : b_.build_vadd(a, b); }
Value* MirLaneEmitter::sub(Value* a, Value* b) { return (lanes_ == 1) ? b_.build_sub(a, b) : b_.build_vsub(a, b); }
Value* MirLaneEmitter::mul(Value* a, Value* b) { return (lanes_ == 1) ? b_.build_mul(a, b) : b_.build_vmul(a, b); }
Value* MirLaneEmitter::fma(Value* a, Value* b, Value* c) { return (lanes_ == 1) ? b_.build_fma(a, b, c) : b_.build_vfma(a, b, c); }
Value* MirLaneEmitter::neg(Value* a) { return (lanes_ == 1) ? b_.build_neg(a) : b_.build_vneg(a); }
Value* MirLaneEmitter::fmin(Value* a, Value* b) { return (lanes_ == 1) ? b_.build_fmin_f32(a, b) : b_.build_vmin(a, b); }
Value* MirLaneEmitter::fmax(Value* a, Value* b) { return (lanes_ == 1) ? b_.build_fmax_f32(a, b) : b_.build_vmax(a, b); }
Value* MirLaneEmitter::sqrt(Value* a) { return (lanes_ == 1) ? b_.build_sqrt_f32(a) : b_.build_vsqrt(a); }
Value* MirLaneEmitter::floor(Value* a) { return (lanes_ == 1) ? b_.build_floor_f32(a) : b_.build_vfloor(a); }
Value* MirLaneEmitter::round(Value* a) { return (lanes_ == 1) ? b_.build_round_f32(a) : b_.build_vround(a); }

Value* MirLaneEmitter::abs(Value* a) {
    return (lanes_ == 1) ? b_.build_fabs_f32(a) : b_.build_vmax(a, b_.build_vneg(a));
}

Value* MirLaneEmitter::rawDiv(Value* a, Value* b) { return (lanes_ == 1) ? b_.build_sdiv(a, b) : b_.build_vdiv(a, b); }

Value* MirLaneEmitter::safeDiv(Value* a, Value* b) {
    return selectGt(abs(b), cF32(1e-20f), rawDiv(a, b), cF32(0.0f));
}

Value* MirLaneEmitter::clamp(Value* x, Value* lo, Value* hi) { return fmin(fmax(x, lo), hi); }
Value* MirLaneEmitter::lerp(Value* a, Value* b, Value* t) { return fma(sub(b, a), t, a); }
Value* MirLaneEmitter::length2(Value* x, Value* y) { return sqrt(fma(y, y, mul(x, x))); }
Value* MirLaneEmitter::length3(Value* x, Value* y, Value* z) { return sqrt(fma(z, z, fma(y, y, mul(x, x)))); }

Value* MirLaneEmitter::selectGt(Value* lhs, Value* rhs, Value* tVal, Value* fVal) {
    if (lanes_ == 1) {
        return b_.build_select(b_.build_sgt(lhs, rhs), tVal, fVal);
    }
    return b_.build_vselect(b_.build_vcmp(Opcode::vcmp_gt, lhs, rhs), tVal, fVal);
}

Value* MirLaneEmitter::fptosiI32(Value* a) {
    return (lanes_ == 1) ? b_.build_fptosi_i32_f32(a) : b_.build_vconvert(i32Type(), a);
}

Value* MirLaneEmitter::sitofpF32(Value* a) {
    return (lanes_ == 1) ? b_.build_sitofp_f32_i32(a) : b_.build_vconvert(f32Type(), a);
}

Value* MirLaneEmitter::bitcastI32(Value* a) {
    return (lanes_ == 1) ? b_.build_bitcast_i32_f32(a) : b_.build_vbitcast(i32Type(), a);
}

Value* MirLaneEmitter::bitcastF32(Value* a) {
    return (lanes_ == 1) ? b_.build_bitcast_f32_i32(a) : b_.build_vbitcast(f32Type(), a);
}

Value* MirLaneEmitter::ixor(Value* a, Value* b) { return (lanes_ == 1) ? b_.build_xor(a, b) : b_.build_vxor(a, b); }
Value* MirLaneEmitter::iand(Value* a, Value* b) { return (lanes_ == 1) ? b_.build_and(a, b) : b_.build_vand(a, b); }
Value* MirLaneEmitter::ishl(Value* a, Value* b) { return (lanes_ == 1) ? b_.build_shl(a, b) : b_.build_vshl(a, b); }
Value* MirLaneEmitter::ilshr(Value* a, Value* b) { return (lanes_ == 1) ? b_.build_lshr(a, b) : b_.build_vlshr(a, b); }

Value* MirLaneEmitter::sinF32(Value* x) {
    if (lanes_ == 1 && gpuMode_) {
        if (b_.current_module() && !b_.current_module()->has_external_symbol("gpu_sin")) {
            b_.current_module()->add_external_symbol("gpu_sin");
        }
        return b_.build_call("gpu_sin", Type::f32(), {x});
    }
    Value* invPi = cF32(0.31830988618379067f);
    Value* k = round(mul(x, invPi));
    Value* negPiHi = cF32(-3.140625f);
    Value* negPiLo = cF32(-0.00096765358979323846f);
    Value* r = fma(k, negPiLo, fma(k, negPiHi, x));
    Value* r2 = mul(r, r);
    Value* r3 = mul(r2, r);
    Value* p = fma(r2, cF32(2.7525562e-6f), cF32(-0.00019840874f));
    p = fma(r2, p, cF32(0.0083333310f));
    p = fma(r2, p, cF32(-0.16666667f));
    p = fma(r3, p, r);
    Value* ki = fptosiI32(k);
    Value* signBit = ishl(iand(ki, cI32(1)), cI32(31));
    return bitcastF32(ixor(bitcastI32(p), signBit));
}

Value* MirLaneEmitter::cosF32(Value* x) {
    if (lanes_ == 1 && gpuMode_) {
        if (b_.current_module() && !b_.current_module()->has_external_symbol("gpu_cos")) {
            b_.current_module()->add_external_symbol("gpu_cos");
        }
        return b_.build_call("gpu_cos", Type::f32(), {x});
    }
    return sinF32(add(x, cF32(1.5707963267948966f)));
}

Value* MirLaneEmitter::loadLaneF32(Value* basePtr, Value* elemIdxI32) {
    Value* byteOff = b_.build_shl(b_.build_sext_i64(elemIdxI32), kb_.const_i64(2));
    Value* ptr = b_.build_add(basePtr, byteOff);
    return (lanes_ == 1) ? b_.build_load(Type::f32(), ptr, 0) : b_.build_vload(f32Type(), ptr, 0);
}

void MirLaneEmitter::storeLaneF32(Value* basePtr, Value* elemIdxI32, Value* val) {
    Value* byteOff = b_.build_shl(b_.build_sext_i64(elemIdxI32), kb_.const_i64(2));
    Value* ptr = b_.build_add(basePtr, byteOff);
    if (lanes_ == 1) {
        b_.build_store(Type::f32(), ptr, 0, val);
    } else {
        b_.build_vstore(f32Type(), ptr, 0, val);
    }
}

Value* MirLaneEmitter::loadScalarBroadcastF32(Value* ptr, int32_t byteOffset) {
    return broadcastF32(b_.build_load(Type::f32(), ptr, byteOffset));
}

Value* buildLaneRampF32(MirLaneEmitter& em) {
    Value* ramp = em.b().build_vzero(em.f32Type());
    for (uint32_t lane = 0; lane < em.lanes(); ++lane) {
        ramp = em.b().build_vinsert_lane(ramp, em.kb().const_f32(static_cast<float>(lane)), lane);
    }
    return ramp;
}

namespace {

Value* emitNoise3dSingleSource(MirLaneEmitter& em, Value* px, Value* py, Value* pz) {
    Value* flrX = em.floor(px);
    Value* flrY = em.floor(py);
    Value* flrZ = em.floor(pz);

    Value* i0 = em.fptosiI32(flrX);
    Value* j0 = em.fptosiI32(flrY);
    Value* k0 = em.fptosiI32(flrZ);
    Value* oneI = em.cI32(1);
    Value* i1 = em.add(i0, oneI);
    Value* j1 = em.add(j0, oneI);
    Value* k1 = em.add(k0, oneI);

    Value* fx = em.sub(px, flrX);
    Value* fy = em.sub(py, flrY);
    Value* fz = em.sub(pz, flrZ);
    Value* c3 = em.cF32(3.0f);
    Value* c2 = em.cF32(2.0f);
    Value* u = em.mul(em.mul(fx, fx), em.sub(c3, em.mul(fx, c2)));
    Value* v = em.mul(em.mul(fy, fy), em.sub(c3, em.mul(fy, c2)));
    Value* w = em.mul(em.mul(fz, fz), em.sub(c3, em.mul(fz, c2)));

    Value* kPrime1 = em.cI32(374761393);
    Value* kPrime2 = em.cI32(668265263);
    Value* kPrime3 = em.cI32(1274126177);
    Value* s13 = em.cI32(13);
    Value* s16 = em.cI32(16);
    Value* mask15 = em.cI32(0x7FFF);
    Value* scale = em.cF32(2.0f / 32767.0f);
    Value* oneF = em.cF32(1.0f);

    auto hashCorner = [&](Value* ix, Value* iy, Value* iz) -> Value* {
        Value* h = em.ixor(ix, em.mul(iy, kPrime1));
        h = em.ixor(h, em.mul(iz, kPrime2));
        h = em.mul(em.ixor(h, em.ilshr(h, s13)), kPrime3);
        h = em.ixor(h, em.ilshr(h, s16));
        Value* hF = em.sitofpF32(em.iand(h, mask15));
        return em.sub(em.mul(hF, scale), oneF);
    };

    Value* x00 = em.lerp(hashCorner(i0, j0, k0), hashCorner(i1, j0, k0), u);
    Value* x10 = em.lerp(hashCorner(i0, j1, k0), hashCorner(i1, j1, k0), u);
    Value* x01 = em.lerp(hashCorner(i0, j0, k1), hashCorner(i1, j0, k1), u);
    Value* x11 = em.lerp(hashCorner(i0, j1, k1), hashCorner(i1, j1, k1), u);
    return em.lerp(em.lerp(x00, x10, v), em.lerp(x01, x11, v), w);
}

Value* emitCustomExprSingleSource(MirLaneEmitter& em, const SdfCustomExpr& expr,
                                  Value* px, Value* py, Value* pz) {
    if (expr.instructions.empty() || expr.resultReg < 0) {
        return em.cF32(1e9f);
    }
    std::vector<Value*> regs(expr.instructions.size(), nullptr);
    Value* zero = em.cF32(0.0f);
    for (size_t i = 0; i < expr.instructions.size(); ++i) {
        const SdfExprInst& inst = expr.instructions[i];
        Value* a = (inst.a >= 0) ? regs[static_cast<size_t>(inst.a)] : nullptr;
        Value* b = (inst.b >= 0) ? regs[static_cast<size_t>(inst.b)] : nullptr;
        Value* c = (inst.c >= 0) ? regs[static_cast<size_t>(inst.c)] : nullptr;
        switch (inst.op) {
            case SdfExprOp::VarX:  regs[i] = px; break;
            case SdfExprOp::VarY:  regs[i] = py; break;
            case SdfExprOp::VarZ:  regs[i] = pz; break;
            case SdfExprOp::Const: regs[i] = em.cF32(inst.imm); break;
            case SdfExprOp::Neg:   regs[i] = em.neg(a); break;
            case SdfExprOp::Add:   regs[i] = em.add(a, b); break;
            case SdfExprOp::Sub:   regs[i] = em.sub(a, b); break;
            case SdfExprOp::Mul:   regs[i] = em.mul(a, b); break;
            case SdfExprOp::Div:   regs[i] = em.safeDiv(a, b); break;
            case SdfExprOp::Sin:   regs[i] = em.sinF32(a); break;
            case SdfExprOp::Cos:   regs[i] = em.cosF32(a); break;
            case SdfExprOp::Abs:   regs[i] = em.abs(a); break;
            case SdfExprOp::Sqrt:  regs[i] = em.sqrt(em.fmax(a, zero)); break;
            case SdfExprOp::Floor: regs[i] = em.floor(a); break;
            case SdfExprOp::Min:   regs[i] = em.fmin(a, b); break;
            case SdfExprOp::Max:   regs[i] = em.fmax(a, b); break;
            case SdfExprOp::Clamp: regs[i] = em.clamp(a, b, c); break;
            case SdfExprOp::Fma:   regs[i] = em.fma(a, b, c); break;
        }
    }
    return regs[static_cast<size_t>(expr.resultReg)];
}

} // namespace

Value* emitSdfGraphSingleSource(MirLaneEmitter& em, const SdfGraph& graph, int nodeId,
                                Value* px, Value* py, Value* pz) {
    if (nodeId < 0 || nodeId >= static_cast<int>(graph.size())) {
        return em.cF32(1e9f);
    }
    const SdfNode& node = graph.node(nodeId);
    Value* zero = em.cF32(0.0f);
    switch (node.op) {
        case SdfOp::Sphere:
            return em.sub(em.length3(px, py, pz), em.cF32(node.f0));
        case SdfOp::Box: {
            Value* qx = em.sub(em.abs(px), em.cF32(node.v0.x));
            Value* qy = em.sub(em.abs(py), em.cF32(node.v0.y));
            Value* qz = em.sub(em.abs(pz), em.cF32(node.v0.z));
            Value* outside = em.length3(em.fmax(qx, zero), em.fmax(qy, zero), em.fmax(qz, zero));
            Value* inside = em.fmin(em.fmax(qx, em.fmax(qy, qz)), zero);
            return em.add(outside, inside);
        }
        case SdfOp::RoundedBox: {
            Value* r = em.cF32(node.f0);
            Value* qx = em.add(em.sub(em.abs(px), em.cF32(node.v0.x)), r);
            Value* qy = em.add(em.sub(em.abs(py), em.cF32(node.v0.y)), r);
            Value* qz = em.add(em.sub(em.abs(pz), em.cF32(node.v0.z)), r);
            Value* outside = em.length3(em.fmax(qx, zero), em.fmax(qy, zero), em.fmax(qz, zero));
            Value* inside = em.fmin(em.fmax(qx, em.fmax(qy, qz)), zero);
            return em.sub(em.add(outside, inside), r);
        }
        case SdfOp::Cylinder: {
            Value* dx = em.sub(em.length2(px, pz), em.cF32(node.f0));
            Value* dy = em.sub(em.abs(py), em.cF32(node.f1));
            Value* outside = em.length2(em.fmax(dx, zero), em.fmax(dy, zero));
            Value* inside = em.fmin(em.fmax(dx, dy), zero);
            return em.add(outside, inside);
        }
        case SdfOp::Capsule: {
            Value* pax = em.sub(px, em.cF32(node.v0.x));
            Value* pay = em.sub(py, em.cF32(node.v0.y));
            Value* paz = em.sub(pz, em.cF32(node.v0.z));
            Value* bax = em.cF32(node.v1.x - node.v0.x);
            Value* bay = em.cF32(node.v1.y - node.v0.y);
            Value* baz = em.cF32(node.v1.z - node.v0.z);
            Value* dotPaBa = em.fma(paz, baz, em.fma(pay, bay, em.mul(pax, bax)));
            Value* dotBaBa = em.fma(baz, baz, em.fma(bay, bay, em.mul(bax, bax)));
            Value* h = em.clamp(em.safeDiv(dotPaBa, dotBaBa), zero, em.cF32(1.0f));
            Value* rx = em.sub(pax, em.mul(bax, h));
            Value* ry = em.sub(pay, em.mul(bay, h));
            Value* rz = em.sub(paz, em.mul(baz, h));
            return em.sub(em.length3(rx, ry, rz), em.cF32(node.f0));
        }
        case SdfOp::Torus: {
            Value* qx = em.sub(em.length2(px, pz), em.cF32(node.f0));
            return em.sub(em.length2(qx, py), em.cF32(node.f1));
        }
        case SdfOp::Plane: {
            Value* dot = em.fma(pz, em.cF32(node.v0.z), em.fma(py, em.cF32(node.v0.y), em.mul(px, em.cF32(node.v0.x))));
            return em.add(dot, em.cF32(node.f0));
        }
        case SdfOp::Union:
            return em.fmin(emitSdfGraphSingleSource(em, graph, node.left, px, py, pz),
                           emitSdfGraphSingleSource(em, graph, node.right, px, py, pz));
        case SdfOp::Intersection:
            return em.fmax(emitSdfGraphSingleSource(em, graph, node.left, px, py, pz),
                           emitSdfGraphSingleSource(em, graph, node.right, px, py, pz));
        case SdfOp::Subtraction:
            return em.fmax(emitSdfGraphSingleSource(em, graph, node.left, px, py, pz),
                           em.neg(emitSdfGraphSingleSource(em, graph, node.right, px, py, pz)));
        case SdfOp::SmoothUnion: {
            Value* d1 = emitSdfGraphSingleSource(em, graph, node.left, px, py, pz);
            Value* d2 = emitSdfGraphSingleSource(em, graph, node.right, px, py, pz);
            Value* k = em.cF32(node.f0);
            Value* half = em.cF32(0.5f);
            Value* one = em.cF32(1.0f);
            Value* h = em.clamp(em.add(half, em.rawDiv(em.mul(half, em.sub(d2, d1)), k)), zero, one);
            Value* mixed = em.lerp(d2, d1, h);
            return em.sub(mixed, em.mul(k, em.mul(h, em.sub(one, h))));
        }
        case SdfOp::SmoothIntersection: {
            Value* d1 = emitSdfGraphSingleSource(em, graph, node.left, px, py, pz);
            Value* d2 = emitSdfGraphSingleSource(em, graph, node.right, px, py, pz);
            Value* k = em.cF32(node.f0);
            Value* half = em.cF32(0.5f);
            Value* one = em.cF32(1.0f);
            Value* h = em.clamp(em.sub(half, em.rawDiv(em.mul(half, em.sub(d2, d1)), k)), zero, one);
            Value* mixed = em.lerp(d2, d1, h);
            return em.add(mixed, em.mul(k, em.mul(h, em.sub(one, h))));
        }
        case SdfOp::SmoothSubtraction: {
            Value* d1 = emitSdfGraphSingleSource(em, graph, node.left, px, py, pz);
            Value* d2 = emitSdfGraphSingleSource(em, graph, node.right, px, py, pz);
            Value* negD2 = em.neg(d2);
            Value* k = em.cF32(node.f0);
            Value* half = em.cF32(0.5f);
            Value* one = em.cF32(1.0f);
            Value* h = em.clamp(em.sub(half, em.rawDiv(em.mul(half, em.add(d2, d1)), k)), zero, one);
            Value* mixed = em.lerp(d1, negD2, h);
            return em.add(mixed, em.mul(k, em.mul(h, em.sub(one, h))));
        }
        case SdfOp::Translate:
            return emitSdfGraphSingleSource(em, graph, node.left,
                                            em.sub(px, em.cF32(node.v0.x)),
                                            em.sub(py, em.cF32(node.v0.y)),
                                            em.sub(pz, em.cF32(node.v0.z)));
        case SdfOp::RotateY: {
            Value* c = em.cF32(std::cos(node.f0));
            Value* s = em.cF32(std::sin(node.f0));
            Value* rx = em.sub(em.mul(px, c), em.mul(pz, s));
            Value* rz = em.add(em.mul(px, s), em.mul(pz, c));
            return emitSdfGraphSingleSource(em, graph, node.left, rx, py, rz);
        }
        case SdfOp::ScaleUniform: {
            Value* s = em.cF32(node.f0);
            Value* dist = emitSdfGraphSingleSource(em, graph, node.left,
                                                   em.rawDiv(px, s), em.rawDiv(py, s), em.rawDiv(pz, s));
            return em.mul(dist, s);
        }
        case SdfOp::Noise3D: {
            Value* n = emitNoise3dSingleSource(em,
                                               em.mul(px, em.cF32(node.v0.x)),
                                               em.mul(py, em.cF32(node.v0.y)),
                                               em.mul(pz, em.cF32(node.v0.z)));
            return em.mul(n, em.cF32(node.f0));
        }
        case SdfOp::Displace:
            return em.add(emitSdfGraphSingleSource(em, graph, node.left, px, py, pz),
                          emitSdfGraphSingleSource(em, graph, node.right, px, py, pz));
        case SdfOp::CustomExpr:
            return emitCustomExprSingleSource(em, graph.customExprData(node.left), px, py, pz);
    }
    return em.cF32(1e9f);
}

} // namespace mir_detail
#endif

} // namespace bromesh
