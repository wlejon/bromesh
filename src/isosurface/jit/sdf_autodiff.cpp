#include "sdf_compiler_internal.h"

#if BROMESH_HAS_BRASS_JIT
#include <brass/mir/gvn.hpp>
#include <brass/mir/scalar_opt.hpp>
#include <brass/mir/smt_guard_elim.hpp>
#include <brass/mir/superopt.hpp>
#endif

#include <algorithm>
#include <cmath>
#include <vector>

namespace bromesh {

#if BROMESH_HAS_BRASS_JIT
using namespace brass;
using namespace brass::codegen;

namespace {

DualF32 dualConst(SdfKernelBuilder& kb, Value* c) {
    Value* z = kb.const_f32(0.0f);
    return {c, z, z, z};
}

DualF32 dualConst(SdfKernelBuilder& kb, float c) {
    return dualConst(kb, kb.const_f32(c));
}

DualF32 dualAdd(SdfKernelBuilder& kb, DualF32 a, DualF32 b) {
    return {kb.add(a.val, b.val), kb.add(a.dx, b.dx), kb.add(a.dy, b.dy), kb.add(a.dz, b.dz)};
}

DualF32 dualSub(SdfKernelBuilder& kb, DualF32 a, DualF32 b) {
    return {kb.sub(a.val, b.val), kb.sub(a.dx, b.dx), kb.sub(a.dy, b.dy), kb.sub(a.dz, b.dz)};
}

DualF32 dualNeg(SdfKernelBuilder& kb, DualF32 a) {
    return {
        kb.builder().build_neg(a.val),
        kb.builder().build_neg(a.dx),
        kb.builder().build_neg(a.dy),
        kb.builder().build_neg(a.dz)
    };
}

DualF32 dualScale(SdfKernelBuilder& kb, DualF32 a, Value* s) {
    return {kb.mul(a.val, s), kb.mul(a.dx, s), kb.mul(a.dy, s), kb.mul(a.dz, s)};
}

DualF32 dualMul(SdfKernelBuilder& kb, DualF32 a, DualF32 b) {
    return {
        kb.mul(a.val, b.val),
        kb.fma(a.val, b.dx, kb.mul(b.val, a.dx)),
        kb.fma(a.val, b.dy, kb.mul(b.val, a.dy)),
        kb.fma(a.val, b.dz, kb.mul(b.val, a.dz))
    };
}

DualF32 dualDiv(SdfKernelBuilder& kb, DualF32 a, DualF32 b) {
    Value* absB = kb.abs_f32(b.val);
    Value* safeB = kb.builder().build_select(
        kb.builder().build_sgt(absB, kb.const_f32(1e-20f)),
        b.val, kb.const_f32(1.0f));
    Value* inv = kb.div(kb.const_f32(1.0f), safeB);
    Value* val = kb.mul(a.val, inv);
    Value* negVal = kb.builder().build_neg(val);
    return {
        val,
        kb.mul(kb.fma(negVal, b.dx, a.dx), inv),
        kb.mul(kb.fma(negVal, b.dy, a.dy), inv),
        kb.mul(kb.fma(negVal, b.dz, a.dz), inv)
    };
}

DualF32 dualAbs(SdfKernelBuilder& kb, DualF32 a) {
    Value* nonNeg = kb.builder().build_sge(a.val, kb.const_f32(0.0f));
    Value* sign = kb.builder().build_select(nonNeg, kb.const_f32(1.0f), kb.const_f32(-1.0f));
    return {kb.abs_f32(a.val), kb.mul(a.dx, sign), kb.mul(a.dy, sign), kb.mul(a.dz, sign)};
}

DualF32 dualSqrt(SdfKernelBuilder& kb, DualF32 a) {
    Value* clamped = kb.max_f32(a.val, kb.const_f32(0.0f));
    Value* r = kb.sqrt_f32(clamped);
    Value* pos = kb.builder().build_sgt(r, kb.const_f32(1e-12f));
    Value* safeR = kb.builder().build_select(pos, r, kb.const_f32(1.0f));
    Value* inv2r = kb.builder().build_select(pos, kb.div(kb.const_f32(0.5f), safeR), kb.const_f32(0.0f));
    return {r, kb.mul(a.dx, inv2r), kb.mul(a.dy, inv2r), kb.mul(a.dz, inv2r)};
}

DualF32 dualMin(SdfKernelBuilder& kb, DualF32 a, DualF32 b) {
    Value* c = kb.builder().build_sle(a.val, b.val);
    return {
        kb.min_f32(a.val, b.val),
        kb.builder().build_select(c, a.dx, b.dx),
        kb.builder().build_select(c, a.dy, b.dy),
        kb.builder().build_select(c, a.dz, b.dz)
    };
}

DualF32 dualMax(SdfKernelBuilder& kb, DualF32 a, DualF32 b) {
    Value* c = kb.builder().build_sge(a.val, b.val);
    return {
        kb.max_f32(a.val, b.val),
        kb.builder().build_select(c, a.dx, b.dx),
        kb.builder().build_select(c, a.dy, b.dy),
        kb.builder().build_select(c, a.dz, b.dz)
    };
}

DualF32 dualClamp(SdfKernelBuilder& kb, DualF32 x, DualF32 lo, DualF32 hi) {
    return dualMax(kb, dualMin(kb, x, hi), lo);
}

DualF32 dualLength2(SdfKernelBuilder& kb, DualF32 x, DualF32 y) {
    Value* len = kb.length2(x.val, y.val);
    Value* pos = kb.builder().build_sgt(len, kb.const_f32(1e-12f));
    Value* safeLen = kb.builder().build_select(pos, len, kb.const_f32(1.0f));
    Value* invLen = kb.builder().build_select(pos, kb.div(kb.const_f32(1.0f), safeLen), kb.const_f32(0.0f));
    Value* ux = kb.mul(x.val, invLen);
    Value* uy = kb.mul(y.val, invLen);
    return {
        len,
        kb.fma(ux, x.dx, kb.mul(uy, y.dx)),
        kb.fma(ux, x.dy, kb.mul(uy, y.dy)),
        kb.fma(ux, x.dz, kb.mul(uy, y.dz))
    };
}

DualF32 dualLength3(SdfKernelBuilder& kb, DualF32 x, DualF32 y, DualF32 z) {
    Value* len = kb.length3(x.val, y.val, z.val);
    Value* pos = kb.builder().build_sgt(len, kb.const_f32(1e-12f));
    Value* safeLen = kb.builder().build_select(pos, len, kb.const_f32(1.0f));
    Value* invLen = kb.builder().build_select(pos, kb.div(kb.const_f32(1.0f), safeLen), kb.const_f32(0.0f));
    Value* ux = kb.mul(x.val, invLen);
    Value* uy = kb.mul(y.val, invLen);
    Value* uz = kb.mul(z.val, invLen);
    return {
        len,
        kb.fma(uz, z.dx, kb.fma(uy, y.dx, kb.mul(ux, x.dx))),
        kb.fma(uz, z.dy, kb.fma(uy, y.dy, kb.mul(ux, x.dy))),
        kb.fma(uz, z.dz, kb.fma(uy, y.dz, kb.mul(ux, x.dz)))
    };
}

DualF32 lowerNoise3dDual(SdfKernelBuilder& kb, DualF32 px, DualF32 py, DualF32 pz) {
    Value* flr_x = kb.floor_f32(px.val);
    Value* flr_y = kb.floor_f32(py.val);
    Value* flr_z = kb.floor_f32(pz.val);

    Value* i0 = kb.builder().build_fptosi_i32_f32(flr_x);
    Value* j0 = kb.builder().build_fptosi_i32_f32(flr_y);
    Value* k0 = kb.builder().build_fptosi_i32_f32(flr_z);

    Value* one_i = kb.const_i32(1);
    Value* i1 = kb.builder().build_add(i0, one_i);
    Value* j1 = kb.builder().build_add(j0, one_i);
    Value* k1 = kb.builder().build_add(k0, one_i);

    Value* fx = kb.sub(px.val, flr_x);
    Value* fy = kb.sub(py.val, flr_y);
    Value* fz = kb.sub(pz.val, flr_z);

    Value* c_1 = kb.const_f32(1.0f);
    Value* c_2 = kb.const_f32(2.0f);
    Value* c_3 = kb.const_f32(3.0f);
    Value* c_6 = kb.const_f32(6.0f);

    Value* u = kb.mul(kb.mul(fx, fx), kb.sub(c_3, kb.mul(fx, c_2)));
    Value* v = kb.mul(kb.mul(fy, fy), kb.sub(c_3, kb.mul(fy, c_2)));
    Value* w = kb.mul(kb.mul(fz, fz), kb.sub(c_3, kb.mul(fz, c_2)));

    // du/dfx = 6 * fx * (1 - fx)
    Value* du_dfx = kb.mul(c_6, kb.mul(fx, kb.sub(c_1, fx)));
    Value* dv_dfy = kb.mul(c_6, kb.mul(fy, kb.sub(c_1, fy)));
    Value* dw_dfz = kb.mul(c_6, kb.mul(fz, kb.sub(c_1, fz)));

    auto hash_corner = [&](Value* ix, Value* iy, Value* iz) -> Value* {
        Value* h = kb.builder().build_xor(ix, kb.builder().build_mul(iy, kb.const_i32(374761393)));
        h = kb.builder().build_xor(h, kb.builder().build_mul(iz, kb.const_i32(668265263)));
        Value* shift13 = kb.builder().build_lshr(h, kb.const_i32(13));
        h = kb.builder().build_mul(kb.builder().build_xor(h, shift13), kb.const_i32(1274126177));
        Value* shift16 = kb.builder().build_lshr(h, kb.const_i32(16));
        h = kb.builder().build_xor(h, shift16);
        Value* masked = kb.builder().build_and(h, kb.const_i32(0x7FFF));
        Value* h_f = kb.builder().build_sitofp_f32_i32(masked);
        return kb.sub(kb.mul(h_f, kb.const_f32(2.0f / 32767.0f)), c_1);
    };

    Value* h000 = hash_corner(i0, j0, k0);
    Value* h100 = hash_corner(i1, j0, k0);
    Value* h010 = hash_corner(i0, j1, k0);
    Value* h110 = hash_corner(i1, j1, k0);
    Value* h001 = hash_corner(i0, j0, k1);
    Value* h101 = hash_corner(i1, j0, k1);
    Value* h011 = hash_corner(i0, j1, k1);
    Value* h111 = hash_corner(i1, j1, k1);

    Value* x00 = kb.lerp_f32(h000, h100, u);
    Value* x10 = kb.lerp_f32(h010, h110, u);
    Value* x01 = kb.lerp_f32(h001, h101, u);
    Value* x11 = kb.lerp_f32(h011, h111, u);

    Value* y0 = kb.lerp_f32(x00, x10, v);
    Value* y1 = kb.lerp_f32(x01, x11, v);
    Value* val = kb.lerp_f32(y0, y1, w);

    // Analytic partials wrt (u, v, w)
    Value* dx00 = kb.sub(h100, h000);
    Value* dx10 = kb.sub(h110, h010);
    Value* dx01 = kb.sub(h101, h001);
    Value* dx11 = kb.sub(h111, h011);
    Value* dn_du = kb.lerp_f32(kb.lerp_f32(dx00, dx10, v), kb.lerp_f32(dx01, dx11, v), w);
    Value* dn_dv = kb.lerp_f32(kb.sub(x10, x00), kb.sub(x11, x01), w);
    Value* dn_dw = kb.sub(y1, y0);

    Value* gx_loc = kb.mul(dn_du, du_dfx);
    Value* gy_loc = kb.mul(dn_dv, dv_dfy);
    Value* gz_loc = kb.mul(dn_dw, dw_dfz);

    return {
        val,
        kb.fma(gz_loc, pz.dx, kb.fma(gy_loc, py.dx, kb.mul(gx_loc, px.dx))),
        kb.fma(gz_loc, pz.dy, kb.fma(gy_loc, py.dy, kb.mul(gx_loc, px.dy))),
        kb.fma(gz_loc, pz.dz, kb.fma(gy_loc, py.dz, kb.mul(gx_loc, px.dz)))
    };
}

} // namespace

DualF32 lowerCustomExprDualToMir(SdfKernelBuilder& kb, const SdfCustomExpr& expr,
                                 DualF32 px, DualF32 py, DualF32 pz) {
    if (expr.instructions.empty() || expr.resultReg < 0) {
        return dualConst(kb, 1e9f);
    }
    std::vector<DualF32> regs(expr.instructions.size());
    for (size_t i = 0; i < expr.instructions.size(); ++i) {
        const SdfExprInst& inst = expr.instructions[i];
        DualF32 a = (inst.a >= 0) ? regs[static_cast<size_t>(inst.a)] : dualConst(kb, 0.0f);
        DualF32 b = (inst.b >= 0) ? regs[static_cast<size_t>(inst.b)] : dualConst(kb, 0.0f);
        DualF32 c = (inst.c >= 0) ? regs[static_cast<size_t>(inst.c)] : dualConst(kb, 0.0f);
        switch (inst.op) {
            case SdfExprOp::VarX:  regs[i] = px; break;
            case SdfExprOp::VarY:  regs[i] = py; break;
            case SdfExprOp::VarZ:  regs[i] = pz; break;
            case SdfExprOp::Const: regs[i] = dualConst(kb, inst.imm); break;
            case SdfExprOp::Neg:   regs[i] = dualNeg(kb, a); break;
            case SdfExprOp::Add:   regs[i] = dualAdd(kb, a, b); break;
            case SdfExprOp::Sub:   regs[i] = dualSub(kb, a, b); break;
            case SdfExprOp::Mul:   regs[i] = dualMul(kb, a, b); break;
            case SdfExprOp::Div:   regs[i] = dualDiv(kb, a, b); break;
            case SdfExprOp::Sin: {
                Value* s = emitSinF32(kb, a.val, false);
                Value* co = emitCosF32(kb, a.val, false);
                regs[i] = {s, kb.mul(co, a.dx), kb.mul(co, a.dy), kb.mul(co, a.dz)};
                break;
            }
            case SdfExprOp::Cos: {
                Value* s = emitSinF32(kb, a.val, false);
                Value* co = emitCosF32(kb, a.val, false);
                Value* ns = kb.builder().build_neg(s);
                regs[i] = {co, kb.mul(ns, a.dx), kb.mul(ns, a.dy), kb.mul(ns, a.dz)};
                break;
            }
            case SdfExprOp::Abs:   regs[i] = dualAbs(kb, a); break;
            case SdfExprOp::Sqrt:  regs[i] = dualSqrt(kb, a); break;
            case SdfExprOp::Floor: regs[i] = dualConst(kb, kb.floor_f32(a.val)); break;
            case SdfExprOp::Min:   regs[i] = dualMin(kb, a, b); break;
            case SdfExprOp::Max:   regs[i] = dualMax(kb, a, b); break;
            case SdfExprOp::Clamp: regs[i] = dualClamp(kb, a, b, c); break;
            case SdfExprOp::Fma:   regs[i] = dualAdd(kb, dualMul(kb, a, b), c); break;
        }
    }
    return regs[static_cast<size_t>(expr.resultReg)];
}

DualF32 lowerDualToMir(SdfKernelBuilder& kb, const SdfGraph& graph, int nodeId,
                       DualF32 px, DualF32 py, DualF32 pz) {
    if (nodeId < 0 || nodeId >= static_cast<int>(graph.size())) {
        return dualConst(kb, 1e9f);
    }
    const SdfNode& node = graph.node(nodeId);
    DualF32 zero = dualConst(kb, 0.0f);

    switch (node.op) {
        case SdfOp::Sphere:
            return dualSub(kb, dualLength3(kb, px, py, pz), dualConst(kb, node.f0));

        case SdfOp::Box: {
            DualF32 qx = dualSub(kb, dualAbs(kb, px), dualConst(kb, node.v0.x));
            DualF32 qy = dualSub(kb, dualAbs(kb, py), dualConst(kb, node.v0.y));
            DualF32 qz = dualSub(kb, dualAbs(kb, pz), dualConst(kb, node.v0.z));
            DualF32 outside = dualLength3(kb, dualMax(kb, qx, zero), dualMax(kb, qy, zero), dualMax(kb, qz, zero));
            DualF32 inside = dualMin(kb, dualMax(kb, qx, dualMax(kb, qy, qz)), zero);
            return dualAdd(kb, outside, inside);
        }

        case SdfOp::RoundedBox: {
            DualF32 r = dualConst(kb, node.f0);
            DualF32 qx = dualAdd(kb, dualSub(kb, dualAbs(kb, px), dualConst(kb, node.v0.x)), r);
            DualF32 qy = dualAdd(kb, dualSub(kb, dualAbs(kb, py), dualConst(kb, node.v0.y)), r);
            DualF32 qz = dualAdd(kb, dualSub(kb, dualAbs(kb, pz), dualConst(kb, node.v0.z)), r);
            DualF32 outside = dualLength3(kb, dualMax(kb, qx, zero), dualMax(kb, qy, zero), dualMax(kb, qz, zero));
            DualF32 inside = dualMin(kb, dualMax(kb, qx, dualMax(kb, qy, qz)), zero);
            return dualSub(kb, dualAdd(kb, outside, inside), r);
        }

        case SdfOp::Cylinder: {
            DualF32 dx = dualSub(kb, dualLength2(kb, px, pz), dualConst(kb, node.f0));
            DualF32 dy = dualSub(kb, dualAbs(kb, py), dualConst(kb, node.f1));
            DualF32 outside = dualLength2(kb, dualMax(kb, dx, zero), dualMax(kb, dy, zero));
            DualF32 inside = dualMin(kb, dualMax(kb, dx, dy), zero);
            return dualAdd(kb, outside, inside);
        }

        case SdfOp::Capsule: {
            DualF32 pax = dualSub(kb, px, dualConst(kb, node.v0.x));
            DualF32 pay = dualSub(kb, py, dualConst(kb, node.v0.y));
            DualF32 paz = dualSub(kb, pz, dualConst(kb, node.v0.z));
            float bax = node.v1.x - node.v0.x;
            float bay = node.v1.y - node.v0.y;
            float baz = node.v1.z - node.v0.z;
            float baLen2 = std::max(bax * bax + bay * bay + baz * baz, 1e-12f);
            DualF32 dotPaBa = dualAdd(kb, dualScale(kb, pax, kb.const_f32(bax)),
                              dualAdd(kb, dualScale(kb, pay, kb.const_f32(bay)),
                                          dualScale(kb, paz, kb.const_f32(baz))));
            DualF32 h = dualClamp(kb, dualScale(kb, dotPaBa, kb.const_f32(1.0f / baLen2)), zero, dualConst(kb, 1.0f));
            DualF32 vx = dualSub(kb, pax, dualScale(kb, h, kb.const_f32(bax)));
            DualF32 vy = dualSub(kb, pay, dualScale(kb, h, kb.const_f32(bay)));
            DualF32 vz = dualSub(kb, paz, dualScale(kb, h, kb.const_f32(baz)));
            return dualSub(kb, dualLength3(kb, vx, vy, vz), dualConst(kb, node.f0));
        }

        case SdfOp::Torus: {
            DualF32 qx = dualSub(kb, dualLength2(kb, px, pz), dualConst(kb, node.f0));
            return dualSub(kb, dualLength2(kb, qx, py), dualConst(kb, node.f1));
        }

        case SdfOp::Plane: {
            DualF32 dot = dualAdd(kb, dualScale(kb, px, kb.const_f32(node.v0.x)),
                          dualAdd(kb, dualScale(kb, py, kb.const_f32(node.v0.y)),
                                      dualScale(kb, pz, kb.const_f32(node.v0.z))));
            return dualAdd(kb, dot, dualConst(kb, node.f0));
        }

        case SdfOp::Union:
            return dualMin(kb, lowerDualToMir(kb, graph, node.left, px, py, pz),
                               lowerDualToMir(kb, graph, node.right, px, py, pz));

        case SdfOp::Intersection:
            return dualMax(kb, lowerDualToMir(kb, graph, node.left, px, py, pz),
                               lowerDualToMir(kb, graph, node.right, px, py, pz));

        case SdfOp::Subtraction:
            return dualMax(kb, lowerDualToMir(kb, graph, node.left, px, py, pz),
                               dualNeg(kb, lowerDualToMir(kb, graph, node.right, px, py, pz)));

        case SdfOp::SmoothUnion: {
            DualF32 d1 = lowerDualToMir(kb, graph, node.left, px, py, pz);
            DualF32 d2 = lowerDualToMir(kb, graph, node.right, px, py, pz);
            Value* k = kb.const_f32(std::max(node.f0, 1e-7f));
            Value* half = kb.const_f32(0.5f);
            Value* one = kb.const_f32(1.0f);
            Value* h = kb.clamp_f32(kb.add(half, kb.div(kb.mul(half, kb.sub(d2.val, d1.val)), k)),
                                    kb.const_f32(0.0f), one);
            Value* omh = kb.sub(one, h);
            Value* val = kb.sub(kb.add(kb.mul(d2.val, omh), kb.mul(d1.val, h)), kb.mul(k, kb.mul(h, omh)));
            // Exact partial derivatives: d(smin)/d(d1) = h, d(smin)/d(d2) = 1 - h
            return {
                val,
                kb.fma(d1.dx, h, kb.mul(d2.dx, omh)),
                kb.fma(d1.dy, h, kb.mul(d2.dy, omh)),
                kb.fma(d1.dz, h, kb.mul(d2.dz, omh))
            };
        }

        case SdfOp::SmoothIntersection: {
            DualF32 d1 = lowerDualToMir(kb, graph, node.left, px, py, pz);
            DualF32 d2 = lowerDualToMir(kb, graph, node.right, px, py, pz);
            Value* k = kb.const_f32(std::max(node.f0, 1e-7f));
            Value* half = kb.const_f32(0.5f);
            Value* one = kb.const_f32(1.0f);
            Value* h = kb.clamp_f32(kb.sub(half, kb.div(kb.mul(half, kb.sub(d2.val, d1.val)), k)),
                                    kb.const_f32(0.0f), one);
            Value* omh = kb.sub(one, h);
            Value* val = kb.add(kb.add(kb.mul(d2.val, omh), kb.mul(d1.val, h)), kb.mul(k, kb.mul(h, omh)));
            return {
                val,
                kb.fma(d1.dx, h, kb.mul(d2.dx, omh)),
                kb.fma(d1.dy, h, kb.mul(d2.dy, omh)),
                kb.fma(d1.dz, h, kb.mul(d2.dz, omh))
            };
        }

        case SdfOp::SmoothSubtraction: {
            DualF32 d1 = lowerDualToMir(kb, graph, node.left, px, py, pz);
            DualF32 d2 = dualNeg(kb, lowerDualToMir(kb, graph, node.right, px, py, pz));
            Value* k = kb.const_f32(std::max(node.f0, 1e-7f));
            Value* half = kb.const_f32(0.5f);
            Value* one = kb.const_f32(1.0f);
            Value* h = kb.clamp_f32(kb.sub(half, kb.div(kb.mul(half, kb.sub(d2.val, d1.val)), k)),
                                    kb.const_f32(0.0f), one);
            Value* omh = kb.sub(one, h);
            Value* val = kb.add(kb.add(kb.mul(d2.val, omh), kb.mul(d1.val, h)), kb.mul(k, kb.mul(h, omh)));
            return {
                val,
                kb.fma(d1.dx, h, kb.mul(d2.dx, omh)),
                kb.fma(d1.dy, h, kb.mul(d2.dy, omh)),
                kb.fma(d1.dz, h, kb.mul(d2.dz, omh))
            };
        }

        case SdfOp::Translate: {
            DualF32 tx = dualSub(kb, px, dualConst(kb, node.v0.x));
            DualF32 ty = dualSub(kb, py, dualConst(kb, node.v0.y));
            DualF32 tz = dualSub(kb, pz, dualConst(kb, node.v0.z));
            return lowerDualToMir(kb, graph, node.left, tx, ty, tz);
        }

        case SdfOp::RotateY: {
            Value* c = kb.const_f32(std::cos(node.f0));
            Value* s = kb.const_f32(std::sin(node.f0));
            DualF32 rx = dualSub(kb, dualScale(kb, px, c), dualScale(kb, pz, s));
            DualF32 rz = dualAdd(kb, dualScale(kb, px, s), dualScale(kb, pz, c));
            return lowerDualToMir(kb, graph, node.left, rx, py, rz);
        }

        case SdfOp::ScaleUniform: {
            float sVal = (std::fabs(node.f0) > 1e-12f) ? node.f0 : 1.0f;
            Value* invS = kb.const_f32(1.0f / sVal);
            DualF32 d = lowerDualToMir(kb, graph, node.left,
                                       dualScale(kb, px, invS),
                                       dualScale(kb, py, invS),
                                       dualScale(kb, pz, invS));
            return dualScale(kb, d, kb.const_f32(sVal));
        }

        case SdfOp::Noise3D: {
            DualF32 nx = dualScale(kb, px, kb.const_f32(node.v0.x));
            DualF32 ny = dualScale(kb, py, kb.const_f32(node.v0.y));
            DualF32 nz = dualScale(kb, pz, kb.const_f32(node.v0.z));
            return dualScale(kb, lowerNoise3dDual(kb, nx, ny, nz), kb.const_f32(node.f0));
        }

        case SdfOp::Displace:
            return dualAdd(kb, lowerDualToMir(kb, graph, node.left, px, py, pz),
                               lowerDualToMir(kb, graph, node.right, px, py, pz));

        case SdfOp::CustomExpr:
            return lowerCustomExprDualToMir(kb, graph.customExprData(node.left), px, py, pz);
    }
    return dualConst(kb, 1e9f);
}

Function* buildSdfGradFunction(Module& mod, const SdfGraph& graph, int rootNode) {
    int root = rootNode >= 0 ? rootNode : graph.root();
    Function* fn = mod.create_function("eval_point_and_grad", Type::f32(), {
        Type::f32(), // float x
        Type::f32(), // float y
        Type::f32(), // float z
        Type::ptr(), // float* out_gx
        Type::ptr(), // float* out_gy
        Type::ptr()  // float* out_gz
    });

    SdfKernelBuilder kb(mod, fn);
    BasicBlock* entry = kb.builder().append_block("entry");
    Value* px = kb.builder().add_block_param(entry, Type::f32());
    Value* py = kb.builder().add_block_param(entry, Type::f32());
    Value* pz = kb.builder().add_block_param(entry, Type::f32());
    Value* outGx = kb.builder().add_block_param(entry, Type::ptr());
    Value* outGy = kb.builder().add_block_param(entry, Type::ptr());
    Value* outGz = kb.builder().add_block_param(entry, Type::ptr());

    kb.position_at_end(entry);
    Value* one = kb.const_f32(1.0f);
    Value* zero = kb.const_f32(0.0f);
    DualF32 dx{px, one, zero, zero};
    DualF32 dy{py, zero, one, zero};
    DualF32 dz{pz, zero, zero, one};

    DualF32 res = lowerDualToMir(kb, graph, root, dx, dy, dz);
    kb.store_f32(outGx, res.dx, 0);
    kb.store_f32(outGy, res.dy, 0);
    kb.store_f32(outGz, res.dz, 0);
    kb.builder().build_ret(res.val);

    dominator_cse(*fn);
    gvn_function(*fn);
    run_smt_guard_elimination(*fn, mod);
    run_mir_superopt(*fn, mod);
    eliminate_dead_code(*fn);

    return fn;
}
#endif

float evalGradScalarNode(const SdfGraph& graph, int nodeId, float px, float py, float pz,
                         float& gx, float& gy, float& gz) {
    if (nodeId < 0 || nodeId >= static_cast<int>(graph.size())) {
        gx = gy = gz = 0.0f;
        return 1e9f;
    }
    const SdfNode& node = graph.node(nodeId);
    if (node.op == SdfOp::CustomExpr) {
        return evalCustomExprGradScalar(graph.customExprData(node.left), px, py, pz, gx, gy, gz);
    }
    const float eps = 1e-4f;
    const float inv2eps = 0.5f / eps;
    float d = evalScalarNode(graph, nodeId, px, py, pz);
    gx = (evalScalarNode(graph, nodeId, px + eps, py, pz) - evalScalarNode(graph, nodeId, px - eps, py, pz)) * inv2eps;
    gy = (evalScalarNode(graph, nodeId, px, py + eps, pz) - evalScalarNode(graph, nodeId, px, py - eps, pz)) * inv2eps;
    gz = (evalScalarNode(graph, nodeId, px, py, pz + eps) - evalScalarNode(graph, nodeId, px, py, pz - eps)) * inv2eps;
    return d;
}

SdfContact CompiledSdfKernel::queryContact(float x, float y, float z) const {
    SdfContact contact;
    float gx = 0.0f, gy = 0.0f, gz = 0.0f;
    if (gradFn) {
        contact.distance = gradFn(x, y, z, &gx, &gy, &gz);
    } else if (pointFn) {
        const float eps = 1e-4f;
        const float inv2eps = 0.5f / eps;
        contact.distance = pointFn(x, y, z);
        gx = (pointFn(x + eps, y, z) - pointFn(x - eps, y, z)) * inv2eps;
        gy = (pointFn(x, y + eps, z) - pointFn(x, y - eps, z)) * inv2eps;
        gz = (pointFn(x, y, z + eps) - pointFn(x, y, z - eps)) * inv2eps;
    }
    contact.gradient = {gx, gy, gz};
    float len = std::sqrt(gx * gx + gy * gy + gz * gz);
    if (len > 1e-7f) {
        float invLen = 1.0f / len;
        contact.normal = {gx * invLen, gy * invLen, gz * invLen};
    } else {
        contact.normal = {0.0f, 1.0f, 0.0f};
    }
    contact.closestSurfacePoint = {
        x - contact.distance * contact.normal.x,
        y - contact.distance * contact.normal.y,
        z - contact.distance * contact.normal.z
    };
    return contact;
}

} // namespace bromesh
