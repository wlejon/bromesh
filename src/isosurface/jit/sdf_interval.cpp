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

namespace {

SdfInterval ivAdd(SdfInterval a, SdfInterval b) { return {a.lo + b.lo, a.hi + b.hi}; }
SdfInterval ivSub(SdfInterval a, SdfInterval b) { return {a.lo - b.hi, a.hi - b.lo}; }
SdfInterval ivNeg(SdfInterval a) { return {-a.hi, -a.lo}; }

SdfInterval ivScale(SdfInterval a, float s) {
    return (s >= 0.0f) ? SdfInterval{a.lo * s, a.hi * s} : SdfInterval{a.hi * s, a.lo * s};
}

SdfInterval ivMul(SdfInterval a, SdfInterval b) {
    float p0 = a.lo * b.lo, p1 = a.lo * b.hi, p2 = a.hi * b.lo, p3 = a.hi * b.hi;
    return {std::min({p0, p1, p2, p3}), std::max({p0, p1, p2, p3})};
}

SdfInterval ivDiv(SdfInterval a, SdfInterval b) {
    if (b.lo <= 1e-12f && b.hi >= -1e-12f) {
        return {-1e6f, 1e6f};
    }
    return ivMul(a, {1.0f / b.hi, 1.0f / b.lo});
}

SdfInterval ivAbs(SdfInterval a) {
    return {std::max(0.0f, std::max(a.lo, -a.hi)), std::max(-a.lo, a.hi)};
}

SdfInterval ivSqr(SdfInterval a) {
    SdfInterval ab = ivAbs(a);
    return {ab.lo * ab.lo, ab.hi * ab.hi};
}

SdfInterval ivSqrt(SdfInterval a) {
    return {std::sqrt(std::max(0.0f, a.lo)), std::sqrt(std::max(0.0f, a.hi))};
}

SdfInterval ivMin(SdfInterval a, SdfInterval b) { return {std::min(a.lo, b.lo), std::min(a.hi, b.hi)}; }
SdfInterval ivMax(SdfInterval a, SdfInterval b) { return {std::max(a.lo, b.lo), std::max(a.hi, b.hi)}; }

SdfInterval ivClamp(SdfInterval x, SdfInterval lo, SdfInterval hi) {
    return ivMax(ivMin(x, hi), lo);
}

SdfInterval ivSin(SdfInterval a) {
    const float twoPi = 6.28318530717958647692f;
    const float halfPi = 1.57079632679489661923f;
    if (a.hi - a.lo >= twoPi) return {-1.0f, 1.0f};
    float s0 = std::sin(a.lo), s1 = std::sin(a.hi);
    float lo = std::min(s0, s1), hi = std::max(s0, s1);
    // Maxima at pi/2 + 2*pi*k -> (a - pi/2) / (2*pi) crosses an integer
    if (std::floor((a.lo - halfPi) / twoPi) < std::floor((a.hi - halfPi) / twoPi)) hi = 1.0f;
    // Minima at -pi/2 + 2*pi*k -> (a + pi/2) / (2*pi) crosses an integer
    if (std::floor((a.lo + halfPi) / twoPi) < std::floor((a.hi + halfPi) / twoPi)) lo = -1.0f;
    return {lo, hi};
}

SdfInterval ivCos(SdfInterval a) {
    const float halfPi = 1.57079632679489661923f;
    return ivSin({a.lo + halfPi, a.hi + halfPi});
}

SdfInterval ivLength2(SdfInterval x, SdfInterval y) {
    return ivSqrt(ivAdd(ivSqr(x), ivSqr(y)));
}

SdfInterval ivLength3(SdfInterval x, SdfInterval y, SdfInterval z) {
    return ivSqrt(ivAdd(ivAdd(ivSqr(x), ivSqr(y)), ivSqr(z)));
}

float scalarBoxSdfFromQ(float qx, float qy, float qz) {
    float ox = std::max(qx, 0.0f), oy = std::max(qy, 0.0f), oz = std::max(qz, 0.0f);
    return std::sqrt(ox * ox + oy * oy + oz * oz) + std::min(std::max({qx, qy, qz}), 0.0f);
}

float scalarCylSdfFromD(float dx, float dy) {
    float ox = std::max(dx, 0.0f), oy = std::max(dy, 0.0f);
    return std::sqrt(ox * ox + oy * oy) + std::min(std::max(dx, dy), 0.0f);
}

float scalarSmoothUnion(float d1, float d2, float k) {
    float sk = std::max(k, 1e-7f);
    float h = std::clamp(0.5f + 0.5f * (d2 - d1) / sk, 0.0f, 1.0f);
    return (d2 + (d1 - d2) * h) - sk * h * (1.0f - h);
}

float scalarSmoothIntersection(float d1, float d2, float k) {
    float sk = std::max(k, 1e-7f);
    float h = std::clamp(0.5f - 0.5f * (d2 - d1) / sk, 0.0f, 1.0f);
    return (d2 + (d1 - d2) * h) + sk * h * (1.0f - h);
}

SdfInterval evalNodeIntervalImpl(const SdfGraph& graph, int nodeId,
                                 SdfInterval bx, SdfInterval by, SdfInterval bz,
                                 int* activeNodeOut) {
    if (activeNodeOut) *activeNodeOut = nodeId;
    if (nodeId < 0 || nodeId >= static_cast<int>(graph.size())) {
        return {1e9f, 1e9f};
    }
    const SdfNode& node = graph.node(nodeId);
    switch (node.op) {
        case SdfOp::Sphere: {
            SdfInterval l = ivLength3(bx, by, bz);
            return {l.lo - node.f0, l.hi - node.f0};
        }
        case SdfOp::Box: {
            SdfInterval qx = ivSub(ivAbs(bx), {node.v0.x, node.v0.x});
            SdfInterval qy = ivSub(ivAbs(by), {node.v0.y, node.v0.y});
            SdfInterval qz = ivSub(ivAbs(bz), {node.v0.z, node.v0.z});
            return {scalarBoxSdfFromQ(qx.lo, qy.lo, qz.lo), scalarBoxSdfFromQ(qx.hi, qy.hi, qz.hi)};
        }
        case SdfOp::RoundedBox: {
            float r = node.f0;
            SdfInterval qx = ivSub(ivAbs(bx), {node.v0.x - r, node.v0.x - r});
            SdfInterval qy = ivSub(ivAbs(by), {node.v0.y - r, node.v0.y - r});
            SdfInterval qz = ivSub(ivAbs(bz), {node.v0.z - r, node.v0.z - r});
            return {scalarBoxSdfFromQ(qx.lo, qy.lo, qz.lo) - r, scalarBoxSdfFromQ(qx.hi, qy.hi, qz.hi) - r};
        }
        case SdfOp::Cylinder: {
            SdfInterval dx = ivSub(ivLength2(bx, bz), {node.f0, node.f0});
            SdfInterval dy = ivSub(ivAbs(by), {node.f1, node.f1});
            return {scalarCylSdfFromD(dx.lo, dy.lo), scalarCylSdfFromD(dx.hi, dy.hi)};
        }
        case SdfOp::Capsule: {
            SdfInterval pax = ivSub(bx, {node.v0.x, node.v0.x});
            SdfInterval pay = ivSub(by, {node.v0.y, node.v0.y});
            SdfInterval paz = ivSub(bz, {node.v0.z, node.v0.z});
            float bax = node.v1.x - node.v0.x, bay = node.v1.y - node.v0.y, baz = node.v1.z - node.v0.z;
            float baLen2 = std::max(bax * bax + bay * bay + baz * baz, 1e-12f);
            SdfInterval dot = ivAdd(ivAdd(ivScale(pax, bax), ivScale(pay, bay)), ivScale(paz, baz));
            SdfInterval h = ivClamp(ivScale(dot, 1.0f / baLen2), {0.0f, 0.0f}, {1.0f, 1.0f});
            SdfInterval vx = ivSub(pax, ivScale(h, bax));
            SdfInterval vy = ivSub(pay, ivScale(h, bay));
            SdfInterval vz = ivSub(paz, ivScale(h, baz));
            SdfInterval l = ivLength3(vx, vy, vz);
            SdfInterval nat{l.lo - node.f0, l.hi - node.f0};
            // Intersect with 1-Lipschitz center bound for extra tightness
            float cx = 0.5f * (bx.lo + bx.hi), cy = 0.5f * (by.lo + by.hi), cz = 0.5f * (bz.lo + bz.hi);
            float ex = 0.5f * (bx.hi - bx.lo), ey = 0.5f * (by.hi - by.lo), ez = 0.5f * (bz.hi - bz.lo);
            float rad = std::sqrt(ex * ex + ey * ey + ez * ez);
            float dc = evalScalarNode(graph, nodeId, cx, cy, cz);
            return {std::max(nat.lo, dc - rad), std::min(nat.hi, dc + rad)};
        }
        case SdfOp::Torus: {
            SdfInterval qx = ivAbs(ivSub(ivLength2(bx, bz), {node.f0, node.f0}));
            SdfInterval qy = ivAbs(by);
            SdfInterval l = ivLength2(qx, qy);
            return {l.lo - node.f1, l.hi - node.f1};
        }
        case SdfOp::Plane: {
            SdfInterval dot = ivAdd(ivAdd(ivScale(bx, node.v0.x), ivScale(by, node.v0.y)), ivScale(bz, node.v0.z));
            return {dot.lo + node.f0, dot.hi + node.f0};
        }
        case SdfOp::Union: {
            int actL = node.left, actR = node.right;
            SdfInterval a = evalNodeIntervalImpl(graph, node.left, bx, by, bz, &actL);
            SdfInterval b = evalNodeIntervalImpl(graph, node.right, bx, by, bz, &actR);
            if (activeNodeOut) {
                if (a.hi < b.lo) *activeNodeOut = actL;
                else if (b.hi < a.lo) *activeNodeOut = actR;
            }
            return ivMin(a, b);
        }
        case SdfOp::Intersection: {
            int actL = node.left, actR = node.right;
            SdfInterval a = evalNodeIntervalImpl(graph, node.left, bx, by, bz, &actL);
            SdfInterval b = evalNodeIntervalImpl(graph, node.right, bx, by, bz, &actR);
            if (activeNodeOut) {
                if (a.lo > b.hi) *activeNodeOut = actL;
                else if (b.lo > a.hi) *activeNodeOut = actR;
            }
            return ivMax(a, b);
        }
        case SdfOp::Subtraction: {
            int actL = node.left, actR = node.right;
            SdfInterval a = evalNodeIntervalImpl(graph, node.left, bx, by, bz, &actL);
            SdfInterval b = ivNeg(evalNodeIntervalImpl(graph, node.right, bx, by, bz, &actR));
            if (activeNodeOut && a.lo > b.hi) {
                *activeNodeOut = actL;
            }
            return ivMax(a, b);
        }
        case SdfOp::SmoothUnion: {
            int actL = node.left, actR = node.right;
            SdfInterval a = evalNodeIntervalImpl(graph, node.left, bx, by, bz, &actL);
            SdfInterval b = evalNodeIntervalImpl(graph, node.right, bx, by, bz, &actR);
            float k = std::max(node.f0, 0.0f);
            if (activeNodeOut) {
                if (b.lo - a.hi >= k) *activeNodeOut = actL;
                else if (a.lo - b.hi >= k) *activeNodeOut = actR;
            }
            return {scalarSmoothUnion(a.lo, b.lo, k), scalarSmoothUnion(a.hi, b.hi, k)};
        }
        case SdfOp::SmoothIntersection: {
            int actL = node.left, actR = node.right;
            SdfInterval a = evalNodeIntervalImpl(graph, node.left, bx, by, bz, &actL);
            SdfInterval b = evalNodeIntervalImpl(graph, node.right, bx, by, bz, &actR);
            float k = std::max(node.f0, 0.0f);
            if (activeNodeOut) {
                if (a.lo - b.hi >= k) *activeNodeOut = actL;
                else if (b.lo - a.hi >= k) *activeNodeOut = actR;
            }
            return {scalarSmoothIntersection(a.lo, b.lo, k), scalarSmoothIntersection(a.hi, b.hi, k)};
        }
        case SdfOp::SmoothSubtraction: {
            int actL = node.left, actR = node.right;
            SdfInterval a = evalNodeIntervalImpl(graph, node.left, bx, by, bz, &actL);
            SdfInterval b = ivNeg(evalNodeIntervalImpl(graph, node.right, bx, by, bz, &actR));
            float k = std::max(node.f0, 0.0f);
            if (activeNodeOut && a.lo - b.hi >= k) {
                *activeNodeOut = actL;
            }
            return {scalarSmoothIntersection(a.lo, b.lo, k), scalarSmoothIntersection(a.hi, b.hi, k)};
        }
        case SdfOp::Translate:
            return evalNodeIntervalImpl(graph, node.left,
                                        {bx.lo - node.v0.x, bx.hi - node.v0.x},
                                        {by.lo - node.v0.y, by.hi - node.v0.y},
                                        {bz.lo - node.v0.z, bz.hi - node.v0.z},
                                        nullptr);
        case SdfOp::RotateY: {
            float c = std::cos(node.f0), s = std::sin(node.f0);
            SdfInterval rx = ivSub(ivScale(bx, c), ivScale(bz, s));
            SdfInterval rz = ivAdd(ivScale(bx, s), ivScale(bz, c));
            return evalNodeIntervalImpl(graph, node.left, rx, by, rz, nullptr);
        }
        case SdfOp::ScaleUniform: {
            float s = (std::fabs(node.f0) > 1e-12f) ? node.f0 : 1.0f;
            float invS = 1.0f / s;
            SdfInterval childIv = evalNodeIntervalImpl(graph, node.left,
                                                       ivScale(bx, invS), ivScale(by, invS), ivScale(bz, invS),
                                                       nullptr);
            return ivScale(childIv, std::fabs(s));
        }
        case SdfOp::Noise3D: {
            float amp = std::fabs(node.f0);
            return {-amp, amp};
        }
        case SdfOp::Displace: {
            SdfInterval d = evalNodeIntervalImpl(graph, node.left, bx, by, bz, nullptr);
            SdfInterval n = evalNodeIntervalImpl(graph, node.right, bx, by, bz, nullptr);
            return ivAdd(d, n);
        }
        case SdfOp::CustomExpr:
            return evalCustomExprInterval(graph.customExprData(node.left), bx, by, bz);
    }
    return {1e9f, 1e9f};
}

int clonePrunedSubtree(const SdfGraph& src, int nodeId,
                       SdfInterval bx, SdfInterval by, SdfInterval bz,
                       SdfGraph& dst) {
    if (nodeId < 0 || nodeId >= static_cast<int>(src.size())) return -1;
    int domNode = nodeId;
    evalNodeIntervalImpl(src, nodeId, bx, by, bz, &domNode);
    if (domNode != nodeId) {
        return clonePrunedSubtree(src, domNode, bx, by, bz, dst);
    }
    const SdfNode& n = src.node(nodeId);
    switch (n.op) {
        case SdfOp::Sphere:       return dst.sphere(n.f0);
        case SdfOp::Box:          return dst.box(n.v0);
        case SdfOp::RoundedBox:   return dst.roundedBox(n.v0, n.f0);
        case SdfOp::Cylinder:     return dst.cylinder(n.f0, n.f1);
        case SdfOp::Capsule:      return dst.capsule(n.v0, n.v1, n.f0);
        case SdfOp::Torus:        return dst.torus(n.f0, n.f1);
        case SdfOp::Plane:        return dst.plane(n.v0, n.f0);
        case SdfOp::CustomExpr:   return dst.customExpr(src.customExprData(n.left));
        case SdfOp::Noise3D:      return dst.noise3D(n.v0, n.f0);
        case SdfOp::Union:
            return dst.opUnion(clonePrunedSubtree(src, n.left, bx, by, bz, dst),
                               clonePrunedSubtree(src, n.right, bx, by, bz, dst));
        case SdfOp::Intersection:
            return dst.opIntersection(clonePrunedSubtree(src, n.left, bx, by, bz, dst),
                                      clonePrunedSubtree(src, n.right, bx, by, bz, dst));
        case SdfOp::Subtraction:
            return dst.opSubtraction(clonePrunedSubtree(src, n.left, bx, by, bz, dst),
                                     clonePrunedSubtree(src, n.right, bx, by, bz, dst));
        case SdfOp::SmoothUnion:
            return dst.opSmoothUnion(clonePrunedSubtree(src, n.left, bx, by, bz, dst),
                                     clonePrunedSubtree(src, n.right, bx, by, bz, dst), n.f0);
        case SdfOp::SmoothIntersection:
            return dst.opSmoothIntersection(clonePrunedSubtree(src, n.left, bx, by, bz, dst),
                                            clonePrunedSubtree(src, n.right, bx, by, bz, dst), n.f0);
        case SdfOp::SmoothSubtraction:
            return dst.opSmoothSubtraction(clonePrunedSubtree(src, n.left, bx, by, bz, dst),
                                           clonePrunedSubtree(src, n.right, bx, by, bz, dst), n.f0);
        case SdfOp::Translate:
            return dst.translate(clonePrunedSubtree(src, n.left,
                                                    {bx.lo - n.v0.x, bx.hi - n.v0.x},
                                                    {by.lo - n.v0.y, by.hi - n.v0.y},
                                                    {bz.lo - n.v0.z, bz.hi - n.v0.z}, dst), n.v0);
        case SdfOp::RotateY: {
            float c = std::cos(n.f0), s = std::sin(n.f0);
            SdfInterval rx = ivSub(ivScale(bx, c), ivScale(bz, s));
            SdfInterval rz = ivAdd(ivScale(bx, s), ivScale(bz, c));
            return dst.rotateY(clonePrunedSubtree(src, n.left, rx, by, rz, dst), n.f0);
        }
        case SdfOp::ScaleUniform: {
            float s = (std::fabs(n.f0) > 1e-12f) ? n.f0 : 1.0f;
            float invS = 1.0f / s;
            return dst.scaleUniform(clonePrunedSubtree(src, n.left,
                                                       ivScale(bx, invS), ivScale(by, invS), ivScale(bz, invS), dst), n.f0);
        }
        case SdfOp::Displace:
            return dst.displace(clonePrunedSubtree(src, n.left, bx, by, bz, dst),
                                clonePrunedSubtree(src, n.right, bx, by, bz, dst));
    }
    return -1;
}

} // namespace

SdfInterval evalCustomExprInterval(const SdfCustomExpr& expr, SdfInterval bx, SdfInterval by, SdfInterval bz) {
    if (expr.instructions.empty() || expr.resultReg < 0) {
        return {1e9f, 1e9f};
    }
    std::vector<SdfInterval> regs(expr.instructions.size(), {0.0f, 0.0f});
    for (size_t i = 0; i < expr.instructions.size(); ++i) {
        const SdfExprInst& inst = expr.instructions[i];
        SdfInterval a = (inst.a >= 0) ? regs[static_cast<size_t>(inst.a)] : SdfInterval{0.0f, 0.0f};
        SdfInterval b = (inst.b >= 0) ? regs[static_cast<size_t>(inst.b)] : SdfInterval{0.0f, 0.0f};
        SdfInterval c = (inst.c >= 0) ? regs[static_cast<size_t>(inst.c)] : SdfInterval{0.0f, 0.0f};
        switch (inst.op) {
            case SdfExprOp::VarX:  regs[i] = bx; break;
            case SdfExprOp::VarY:  regs[i] = by; break;
            case SdfExprOp::VarZ:  regs[i] = bz; break;
            case SdfExprOp::Const: regs[i] = {inst.imm, inst.imm}; break;
            case SdfExprOp::Neg:   regs[i] = ivNeg(a); break;
            case SdfExprOp::Add:   regs[i] = ivAdd(a, b); break;
            case SdfExprOp::Sub:   regs[i] = ivSub(a, b); break;
            case SdfExprOp::Mul:   regs[i] = ivMul(a, b); break;
            case SdfExprOp::Div:   regs[i] = ivDiv(a, b); break;
            case SdfExprOp::Sin:   regs[i] = ivSin(a); break;
            case SdfExprOp::Cos:   regs[i] = ivCos(a); break;
            case SdfExprOp::Abs:   regs[i] = ivAbs(a); break;
            case SdfExprOp::Sqrt:  regs[i] = ivSqrt(a); break;
            case SdfExprOp::Floor: regs[i] = {std::floor(a.lo), std::floor(a.hi)}; break;
            case SdfExprOp::Min:   regs[i] = ivMin(a, b); break;
            case SdfExprOp::Max:   regs[i] = ivMax(a, b); break;
            case SdfExprOp::Clamp: regs[i] = ivClamp(a, b, c); break;
            case SdfExprOp::Fma:   regs[i] = ivAdd(ivMul(a, b), c); break;
        }
    }
    return regs[static_cast<size_t>(expr.resultReg)];
}

void evalBlockBounds(const SdfGraph& graph, int nodeId, const bromath::AABB3& blockBounds,
                     float& dMin, float& dMax, int* activeNodeOut) {
    int root = nodeId >= 0 ? nodeId : graph.root();
    SdfInterval bx{std::min(blockBounds.min.x, blockBounds.max.x), std::max(blockBounds.min.x, blockBounds.max.x)};
    SdfInterval by{std::min(blockBounds.min.y, blockBounds.max.y), std::max(blockBounds.min.y, blockBounds.max.y)};
    SdfInterval bz{std::min(blockBounds.min.z, blockBounds.max.z), std::max(blockBounds.min.z, blockBounds.max.z)};
    SdfInterval iv = evalNodeIntervalImpl(graph, root, bx, by, bz, activeNodeOut);
    dMin = iv.lo;
    dMax = iv.hi;
}

SdfGraph pruneCsgGraphForBlock(const SdfGraph& graph, const bromath::AABB3& blockBounds, int rootNode) {
    int root = rootNode >= 0 ? rootNode : graph.root();
    SdfInterval bx{std::min(blockBounds.min.x, blockBounds.max.x), std::max(blockBounds.min.x, blockBounds.max.x)};
    SdfInterval by{std::min(blockBounds.min.y, blockBounds.max.y), std::max(blockBounds.min.y, blockBounds.max.y)};
    SdfInterval bz{std::min(blockBounds.min.z, blockBounds.max.z), std::max(blockBounds.min.z, blockBounds.max.z)};
    SdfGraph pruned;
    int newRoot = clonePrunedSubtree(graph, root, bx, by, bz, pruned);
    if (newRoot >= 0) pruned.setRoot(newRoot);
    return pruned;
}

void CompiledSdfKernel::evalBlockBounds(const bromath::AABB3& blockBounds, float& dMin, float& dMax) const {
    if (blockBoundsFn) {
        blockBoundsFn(blockBounds.min.x, blockBounds.min.y, blockBounds.min.z,
                      blockBounds.max.x, blockBounds.max.y, blockBounds.max.z,
                      &dMin, &dMax);
    } else {
        dMin = -1e9f;
        dMax = 1e9f;
    }
}

#if BROMESH_HAS_BRASS_JIT
using namespace brass;
using namespace brass::codegen;

namespace {

MirInterval mirIvConst(SdfKernelBuilder& kb, float v) {
    Value* c = kb.const_f32(v);
    return {c, c};
}

MirInterval mirIvAdd(SdfKernelBuilder& kb, MirInterval a, MirInterval b) {
    return {kb.add(a.lo, b.lo), kb.add(a.hi, b.hi)};
}

MirInterval mirIvSub(SdfKernelBuilder& kb, MirInterval a, MirInterval b) {
    return {kb.sub(a.lo, b.hi), kb.sub(a.hi, b.lo)};
}

MirInterval mirIvNeg(SdfKernelBuilder& kb, MirInterval a) {
    return {kb.builder().build_neg(a.hi), kb.builder().build_neg(a.lo)};
}

MirInterval mirIvScale(SdfKernelBuilder& kb, MirInterval a, float s) {
    Value* sv = kb.const_f32(s);
    if (s >= 0.0f) return {kb.mul(a.lo, sv), kb.mul(a.hi, sv)};
    return {kb.mul(a.hi, sv), kb.mul(a.lo, sv)};
}

MirInterval mirIvAbs(SdfKernelBuilder& kb, MirInterval a) {
    Value* negLo = kb.builder().build_neg(a.lo);
    Value* negHi = kb.builder().build_neg(a.hi);
    Value* lo = kb.max_f32(kb.const_f32(0.0f), kb.max_f32(a.lo, negHi));
    Value* hi = kb.max_f32(negLo, a.hi);
    return {lo, hi};
}

MirInterval mirIvSqr(SdfKernelBuilder& kb, MirInterval a) {
    MirInterval ab = mirIvAbs(kb, a);
    return {kb.mul(ab.lo, ab.lo), kb.mul(ab.hi, ab.hi)};
}

MirInterval mirIvSqrt(SdfKernelBuilder& kb, MirInterval a) {
    Value* z = kb.const_f32(0.0f);
    return {kb.sqrt_f32(kb.max_f32(a.lo, z)), kb.sqrt_f32(kb.max_f32(a.hi, z))};
}

MirInterval mirIvMul(SdfKernelBuilder& kb, MirInterval a, MirInterval b) {
    Value* p0 = kb.mul(a.lo, b.lo);
    Value* p1 = kb.mul(a.lo, b.hi);
    Value* p2 = kb.mul(a.hi, b.lo);
    Value* p3 = kb.mul(a.hi, b.hi);
    return {
        kb.min_f32(kb.min_f32(p0, p1), kb.min_f32(p2, p3)),
        kb.max_f32(kb.max_f32(p0, p1), kb.max_f32(p2, p3))
    };
}

MirInterval mirIvLength2(SdfKernelBuilder& kb, MirInterval x, MirInterval y) {
    return mirIvSqrt(kb, mirIvAdd(kb, mirIvSqr(kb, x), mirIvSqr(kb, y)));
}

MirInterval mirIvLength3(SdfKernelBuilder& kb, MirInterval x, MirInterval y, MirInterval z) {
    return mirIvSqrt(kb, mirIvAdd(kb, mirIvAdd(kb, mirIvSqr(kb, x), mirIvSqr(kb, y)), mirIvSqr(kb, z)));
}

Value* mirBoxFromQ(SdfKernelBuilder& kb, Value* qx, Value* qy, Value* qz) {
    Value* z = kb.const_f32(0.0f);
    Value* outside = kb.length3(kb.max_f32(qx, z), kb.max_f32(qy, z), kb.max_f32(qz, z));
    Value* inside = kb.min_f32(kb.max_f32(qx, kb.max_f32(qy, qz)), z);
    return kb.add(outside, inside);
}

Value* mirCylFromD(SdfKernelBuilder& kb, Value* dx, Value* dy) {
    Value* z = kb.const_f32(0.0f);
    Value* outside = kb.length2(kb.max_f32(dx, z), kb.max_f32(dy, z));
    Value* inside = kb.min_f32(kb.max_f32(dx, dy), z);
    return kb.add(outside, inside);
}

MirInterval mirIvSin(SdfKernelBuilder& kb, MirInterval a) {
    Value* span = kb.sub(a.hi, a.lo);
    Value* wide = kb.builder().build_sge(span, kb.const_f32(6.283185307179586f));
    Value* s0 = emitSinF32(kb, a.lo, false);
    Value* s1 = emitSinF32(kb, a.hi, false);
    Value* lo = kb.min_f32(s0, s1);
    Value* hi = kb.max_f32(s0, s1);
    Value* invTwoPi = kb.const_f32(0.15915494309189535f);
    Value* halfPi = kb.const_f32(1.5707963267948966f);

    Value* maxCross = kb.builder().build_slt(
        kb.floor_f32(kb.mul(kb.sub(a.lo, halfPi), invTwoPi)),
        kb.floor_f32(kb.mul(kb.sub(a.hi, halfPi), invTwoPi)));
    Value* minCross = kb.builder().build_slt(
        kb.floor_f32(kb.mul(kb.add(a.lo, halfPi), invTwoPi)),
        kb.floor_f32(kb.mul(kb.add(a.hi, halfPi), invTwoPi)));

    hi = kb.builder().build_select(kb.builder().build_or(wide, maxCross), kb.const_f32(1.0f), hi);
    lo = kb.builder().build_select(kb.builder().build_or(wide, minCross), kb.const_f32(-1.0f), lo);
    return {lo, hi};
}

} // namespace

MirInterval lowerCustomExprIntervalToMir(SdfKernelBuilder& kb, const SdfCustomExpr& expr,
                                         MirInterval bx, MirInterval by, MirInterval bz) {
    if (expr.instructions.empty() || expr.resultReg < 0) {
        return mirIvConst(kb, 1e9f);
    }
    std::vector<MirInterval> regs(expr.instructions.size());
    for (size_t i = 0; i < expr.instructions.size(); ++i) {
        const SdfExprInst& inst = expr.instructions[i];
        MirInterval a = (inst.a >= 0) ? regs[static_cast<size_t>(inst.a)] : mirIvConst(kb, 0.0f);
        MirInterval b = (inst.b >= 0) ? regs[static_cast<size_t>(inst.b)] : mirIvConst(kb, 0.0f);
        MirInterval c = (inst.c >= 0) ? regs[static_cast<size_t>(inst.c)] : mirIvConst(kb, 0.0f);
        switch (inst.op) {
            case SdfExprOp::VarX:  regs[i] = bx; break;
            case SdfExprOp::VarY:  regs[i] = by; break;
            case SdfExprOp::VarZ:  regs[i] = bz; break;
            case SdfExprOp::Const: regs[i] = mirIvConst(kb, inst.imm); break;
            case SdfExprOp::Neg:   regs[i] = mirIvNeg(kb, a); break;
            case SdfExprOp::Add:   regs[i] = mirIvAdd(kb, a, b); break;
            case SdfExprOp::Sub:   regs[i] = mirIvSub(kb, a, b); break;
            case SdfExprOp::Mul:   regs[i] = mirIvMul(kb, a, b); break;
            case SdfExprOp::Div: {
                Value* invLo = kb.div(kb.const_f32(1.0f), b.hi);
                Value* invHi = kb.div(kb.const_f32(1.0f), b.lo);
                regs[i] = mirIvMul(kb, a, {kb.min_f32(invLo, invHi), kb.max_f32(invLo, invHi)});
                break;
            }
            case SdfExprOp::Sin:   regs[i] = mirIvSin(kb, a); break;
            case SdfExprOp::Cos: {
                Value* hp = kb.const_f32(1.5707963267948966f);
                regs[i] = mirIvSin(kb, {kb.add(a.lo, hp), kb.add(a.hi, hp)});
                break;
            }
            case SdfExprOp::Abs:   regs[i] = mirIvAbs(kb, a); break;
            case SdfExprOp::Sqrt:  regs[i] = mirIvSqrt(kb, a); break;
            case SdfExprOp::Floor: regs[i] = {kb.floor_f32(a.lo), kb.floor_f32(a.hi)}; break;
            case SdfExprOp::Min:   regs[i] = {kb.min_f32(a.lo, b.lo), kb.min_f32(a.hi, b.hi)}; break;
            case SdfExprOp::Max:   regs[i] = {kb.max_f32(a.lo, b.lo), kb.max_f32(a.hi, b.hi)}; break;
            case SdfExprOp::Clamp:
                regs[i] = {kb.clamp_f32(a.lo, b.lo, c.hi), kb.clamp_f32(a.hi, b.lo, c.hi)};
                break;
            case SdfExprOp::Fma:   regs[i] = mirIvAdd(kb, mirIvMul(kb, a, b), c); break;
        }
    }
    return regs[static_cast<size_t>(expr.resultReg)];
}

MirInterval lowerIntervalToMir(SdfKernelBuilder& kb, const SdfGraph& graph, int nodeId,
                               MirInterval bx, MirInterval by, MirInterval bz) {
    if (nodeId < 0 || nodeId >= static_cast<int>(graph.size())) {
        return mirIvConst(kb, 1e9f);
    }
    const SdfNode& node = graph.node(nodeId);
    switch (node.op) {
        case SdfOp::Sphere: {
            MirInterval l = mirIvLength3(kb, bx, by, bz);
            Value* r = kb.const_f32(node.f0);
            return {kb.sub(l.lo, r), kb.sub(l.hi, r)};
        }
        case SdfOp::Box: {
            MirInterval qx = mirIvSub(kb, mirIvAbs(kb, bx), mirIvConst(kb, node.v0.x));
            MirInterval qy = mirIvSub(kb, mirIvAbs(kb, by), mirIvConst(kb, node.v0.y));
            MirInterval qz = mirIvSub(kb, mirIvAbs(kb, bz), mirIvConst(kb, node.v0.z));
            return {mirBoxFromQ(kb, qx.lo, qy.lo, qz.lo), mirBoxFromQ(kb, qx.hi, qy.hi, qz.hi)};
        }
        case SdfOp::RoundedBox: {
            float r = node.f0;
            MirInterval qx = mirIvSub(kb, mirIvAbs(kb, bx), mirIvConst(kb, node.v0.x - r));
            MirInterval qy = mirIvSub(kb, mirIvAbs(kb, by), mirIvConst(kb, node.v0.y - r));
            MirInterval qz = mirIvSub(kb, mirIvAbs(kb, bz), mirIvConst(kb, node.v0.z - r));
            Value* rv = kb.const_f32(r);
            return {kb.sub(mirBoxFromQ(kb, qx.lo, qy.lo, qz.lo), rv),
                    kb.sub(mirBoxFromQ(kb, qx.hi, qy.hi, qz.hi), rv)};
        }
        case SdfOp::Cylinder: {
            MirInterval dx = mirIvSub(kb, mirIvLength2(kb, bx, bz), mirIvConst(kb, node.f0));
            MirInterval dy = mirIvSub(kb, mirIvAbs(kb, by), mirIvConst(kb, node.f1));
            return {mirCylFromD(kb, dx.lo, dy.lo), mirCylFromD(kb, dx.hi, dy.hi)};
        }
        case SdfOp::Capsule: {
            MirInterval pax = mirIvSub(kb, bx, mirIvConst(kb, node.v0.x));
            MirInterval pay = mirIvSub(kb, by, mirIvConst(kb, node.v0.y));
            MirInterval paz = mirIvSub(kb, bz, mirIvConst(kb, node.v0.z));
            float bax = node.v1.x - node.v0.x, bay = node.v1.y - node.v0.y, baz = node.v1.z - node.v0.z;
            float invLen2 = 1.0f / std::max(bax * bax + bay * bay + baz * baz, 1e-12f);
            MirInterval dot = mirIvAdd(kb, mirIvAdd(kb, mirIvScale(kb, pax, bax), mirIvScale(kb, pay, bay)),
                                       mirIvScale(kb, paz, baz));
            MirInterval hScaled = mirIvScale(kb, dot, invLen2);
            Value* z = kb.const_f32(0.0f), *one = kb.const_f32(1.0f);
            MirInterval h{kb.clamp_f32(hScaled.lo, z, one), kb.clamp_f32(hScaled.hi, z, one)};
            MirInterval vx = mirIvSub(kb, pax, mirIvScale(kb, h, bax));
            MirInterval vy = mirIvSub(kb, pay, mirIvScale(kb, h, bay));
            MirInterval vz = mirIvSub(kb, paz, mirIvScale(kb, h, baz));
            MirInterval l = mirIvLength3(kb, vx, vy, vz);
            Value* r = kb.const_f32(node.f0);
            return {kb.sub(l.lo, r), kb.sub(l.hi, r)};
        }
        case SdfOp::Torus: {
            MirInterval qx = mirIvAbs(kb, mirIvSub(kb, mirIvLength2(kb, bx, bz), mirIvConst(kb, node.f0)));
            MirInterval qy = mirIvAbs(kb, by);
            MirInterval l = mirIvLength2(kb, qx, qy);
            Value* minR = kb.const_f32(node.f1);
            return {kb.sub(l.lo, minR), kb.sub(l.hi, minR)};
        }
        case SdfOp::Plane: {
            MirInterval dot = mirIvAdd(kb, mirIvAdd(kb, mirIvScale(kb, bx, node.v0.x), mirIvScale(kb, by, node.v0.y)),
                                       mirIvScale(kb, bz, node.v0.z));
            Value* d = kb.const_f32(node.f0);
            return {kb.add(dot.lo, d), kb.add(dot.hi, d)};
        }
        case SdfOp::Union: {
            MirInterval a = lowerIntervalToMir(kb, graph, node.left, bx, by, bz);
            MirInterval b = lowerIntervalToMir(kb, graph, node.right, bx, by, bz);
            return {kb.min_f32(a.lo, b.lo), kb.min_f32(a.hi, b.hi)};
        }
        case SdfOp::Intersection: {
            MirInterval a = lowerIntervalToMir(kb, graph, node.left, bx, by, bz);
            MirInterval b = lowerIntervalToMir(kb, graph, node.right, bx, by, bz);
            return {kb.max_f32(a.lo, b.lo), kb.max_f32(a.hi, b.hi)};
        }
        case SdfOp::Subtraction: {
            MirInterval a = lowerIntervalToMir(kb, graph, node.left, bx, by, bz);
            MirInterval b = mirIvNeg(kb, lowerIntervalToMir(kb, graph, node.right, bx, by, bz));
            return {kb.max_f32(a.lo, b.lo), kb.max_f32(a.hi, b.hi)};
        }
        case SdfOp::SmoothUnion: {
            MirInterval a = lowerIntervalToMir(kb, graph, node.left, bx, by, bz);
            MirInterval b = lowerIntervalToMir(kb, graph, node.right, bx, by, bz);
            Value* k = kb.const_f32(std::max(node.f0, 1e-7f));
            return {kb.op_smooth_union(a.lo, b.lo, k), kb.op_smooth_union(a.hi, b.hi, k)};
        }
        case SdfOp::SmoothIntersection: {
            MirInterval a = lowerIntervalToMir(kb, graph, node.left, bx, by, bz);
            MirInterval b = lowerIntervalToMir(kb, graph, node.right, bx, by, bz);
            Value* k = kb.const_f32(std::max(node.f0, 1e-7f));
            return {kb.op_smooth_intersection(a.lo, b.lo, k), kb.op_smooth_intersection(a.hi, b.hi, k)};
        }
        case SdfOp::SmoothSubtraction: {
            MirInterval a = lowerIntervalToMir(kb, graph, node.left, bx, by, bz);
            MirInterval b = mirIvNeg(kb, lowerIntervalToMir(kb, graph, node.right, bx, by, bz));
            Value* k = kb.const_f32(std::max(node.f0, 1e-7f));
            return {kb.op_smooth_intersection(a.lo, b.lo, k), kb.op_smooth_intersection(a.hi, b.hi, k)};
        }
        case SdfOp::Translate:
            return lowerIntervalToMir(kb, graph, node.left,
                                      mirIvSub(kb, bx, mirIvConst(kb, node.v0.x)),
                                      mirIvSub(kb, by, mirIvConst(kb, node.v0.y)),
                                      mirIvSub(kb, bz, mirIvConst(kb, node.v0.z)));
        case SdfOp::RotateY: {
            float c = std::cos(node.f0), s = std::sin(node.f0);
            MirInterval rx = mirIvSub(kb, mirIvScale(kb, bx, c), mirIvScale(kb, bz, s));
            MirInterval rz = mirIvAdd(kb, mirIvScale(kb, bx, s), mirIvScale(kb, bz, c));
            return lowerIntervalToMir(kb, graph, node.left, rx, by, rz);
        }
        case SdfOp::ScaleUniform: {
            float s = (std::fabs(node.f0) > 1e-12f) ? node.f0 : 1.0f;
            float invS = 1.0f / s;
            MirInterval c = lowerIntervalToMir(kb, graph, node.left,
                                               mirIvScale(kb, bx, invS),
                                               mirIvScale(kb, by, invS),
                                               mirIvScale(kb, bz, invS));
            return mirIvScale(kb, c, std::fabs(s));
        }
        case SdfOp::Noise3D: {
            float amp = std::fabs(node.f0);
            return {kb.const_f32(-amp), kb.const_f32(amp)};
        }
        case SdfOp::Displace:
            return mirIvAdd(kb, lowerIntervalToMir(kb, graph, node.left, bx, by, bz),
                                lowerIntervalToMir(kb, graph, node.right, bx, by, bz));
        case SdfOp::CustomExpr:
            return lowerCustomExprIntervalToMir(kb, graph.customExprData(node.left), bx, by, bz);
    }
    return mirIvConst(kb, 1e9f);
}

Function* buildSdfBlockBoundsFunction(Module& mod, const SdfGraph& graph, int rootNode) {
    int root = rootNode >= 0 ? rootNode : graph.root();
    Function* fn = mod.create_function("eval_block_bounds", Type::void_type(), {
        Type::f32(), Type::f32(), Type::f32(), // minX, minY, minZ
        Type::f32(), Type::f32(), Type::f32(), // maxX, maxY, maxZ
        Type::ptr(), Type::ptr()               // float* out_dmin, float* out_dmax
    });

    SdfKernelBuilder kb(mod, fn);
    BasicBlock* entry = kb.builder().append_block("entry");
    Value* minX = kb.builder().add_block_param(entry, Type::f32());
    Value* minY = kb.builder().add_block_param(entry, Type::f32());
    Value* minZ = kb.builder().add_block_param(entry, Type::f32());
    Value* maxX = kb.builder().add_block_param(entry, Type::f32());
    Value* maxY = kb.builder().add_block_param(entry, Type::f32());
    Value* maxZ = kb.builder().add_block_param(entry, Type::f32());
    Value* outMin = kb.builder().add_block_param(entry, Type::ptr());
    Value* outMax = kb.builder().add_block_param(entry, Type::ptr());

    kb.position_at_end(entry);
    MirInterval res = lowerIntervalToMir(kb, graph, root, {minX, maxX}, {minY, maxY}, {minZ, maxZ});
    kb.store_f32(outMin, res.lo, 0);
    kb.store_f32(outMax, res.hi, 0);
    kb.builder().build_ret_void();

    dominator_cse(*fn);
    gvn_function(*fn);
    run_smt_guard_elimination(*fn, mod);
    run_mir_superopt(*fn, mod);
    eliminate_dead_code(*fn);

    return fn;
}
#endif

} // namespace bromesh
