#include "sdf_compiler_internal.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

namespace bromesh {

namespace {

struct ExprParser {
    std::string_view src;
    size_t pos = 0;
    std::vector<SdfExprInst> insts;
    int regX = -1;
    int regY = -1;
    int regZ = -1;

    explicit ExprParser(std::string_view s) : src(s) {
        regX = emit({SdfExprOp::VarX, -1, -1, -1, 0.0f});
        regY = emit({SdfExprOp::VarY, -1, -1, -1, 0.0f});
        regZ = emit({SdfExprOp::VarZ, -1, -1, -1, 0.0f});
    }

    int emit(SdfExprInst inst) {
        int id = static_cast<int>(insts.size());
        insts.push_back(inst);
        return id;
    }

    int emitConst(float v) {
        return emit({SdfExprOp::Const, -1, -1, -1, v});
    }

    void skipWs() {
        while (pos < src.size() && std::isspace(static_cast<unsigned char>(src[pos]))) {
            ++pos;
        }
    }

    bool match(char c) {
        skipWs();
        if (pos < src.size() && src[pos] == c) {
            ++pos;
            return true;
        }
        return false;
    }

    bool matchKeyword(std::string_view kw) {
        skipWs();
        if (pos + kw.size() <= src.size() && src.substr(pos, kw.size()) == kw) {
            size_t next = pos + kw.size();
            if (next == src.size() || (!std::isalnum(static_cast<unsigned char>(src[next])) && src[next] != '_')) {
                pos = next;
                return true;
            }
        }
        return false;
    }

    // Optional JS function / arrow wrapper stripping:
    // e.g. "(x, y, z) => ..." or "function(x, y, z) { return ...; }"
    void stripJsFunctionHeader() {
        skipWs();
        size_t arrowPos = src.find("=>");
        if (arrowPos != std::string_view::npos) {
            pos = arrowPos + 2;
        } else if (matchKeyword("function")) {
            skipWs();
            while (pos < src.size() && src[pos] != '(') ++pos;
            if (pos < src.size() && src[pos] == '(') {
                while (pos < src.size() && src[pos] != ')') ++pos;
                if (pos < src.size()) ++pos;
            }
        }
        skipWs();
        if (match('{')) {
            matchKeyword("return");
        }
    }

    int parseExpression() {
        int lhs = parseTerm();
        for (;;) {
            skipWs();
            if (match('+')) {
                int rhs = parseTerm();
                lhs = emit({SdfExprOp::Add, lhs, rhs, -1, 0.0f});
            } else if (match('-')) {
                int rhs = parseTerm();
                lhs = emit({SdfExprOp::Sub, lhs, rhs, -1, 0.0f});
            } else {
                break;
            }
        }
        return lhs;
    }

    int parseTerm() {
        int lhs = parsePower();
        for (;;) {
            skipWs();
            if (match('*')) {
                int rhs = parsePower();
                lhs = emit({SdfExprOp::Mul, lhs, rhs, -1, 0.0f});
            } else if (match('/')) {
                int rhs = parsePower();
                lhs = emit({SdfExprOp::Div, lhs, rhs, -1, 0.0f});
            } else {
                break;
            }
        }
        return lhs;
    }

    int parsePower() {
        int base = parseUnary();
        skipWs();
        if (pos + 1 < src.size() && src[pos] == '*' && src[pos + 1] == '*') {
            pos += 2;
            return parseExponentTail(base);
        }
        if (match('^')) {
            return parseExponentTail(base);
        }
        return base;
    }

    int parseExponentTail(int base) {
        skipWs();
        int expNode = parseUnary();
        const auto& expInst = insts[static_cast<size_t>(expNode)];
        if (expInst.op == SdfExprOp::Const) {
            float e = expInst.imm;
            if (e == 2.0f) return emit({SdfExprOp::Mul, base, base, -1, 0.0f});
            if (e == 3.0f) {
                int sq = emit({SdfExprOp::Mul, base, base, -1, 0.0f});
                return emit({SdfExprOp::Mul, sq, base, -1, 0.0f});
            }
            if (e == 4.0f) {
                int sq = emit({SdfExprOp::Mul, base, base, -1, 0.0f});
                return emit({SdfExprOp::Mul, sq, sq, -1, 0.0f});
            }
            if (e == 0.5f) return emit({SdfExprOp::Sqrt, base, -1, -1, 0.0f});
            if (e == 0.0f) return emitConst(1.0f);
            if (e == 1.0f) return base;
        }
        return emit({SdfExprOp::Mul, base, base, -1, 0.0f});
    }

    int parseUnary() {
        skipWs();
        if (match('+')) return parseUnary();
        if (match('-')) {
            int v = parseUnary();
            return emit({SdfExprOp::Neg, v, -1, -1, 0.0f});
        }
        return parsePrimary();
    }

    int parsePrimary() {
        skipWs();
        if (match('(')) {
            int v = parseExpression();
            if (!match(')')) {
                throw std::invalid_argument("Custom SDF expression: expected closing ')'");
            }
            return v;
        }

        if (pos < src.size() && (std::isdigit(static_cast<unsigned char>(src[pos])) || src[pos] == '.')) {
            const char* start = src.data() + pos;
            char* endPtr = nullptr;
            float val = std::strtof(start, &endPtr);
            if (endPtr == start) {
                throw std::invalid_argument("Custom SDF expression: invalid numeric literal");
            }
            pos += static_cast<size_t>(endPtr - start);
            if (pos < src.size() && (src[pos] == 'f' || src[pos] == 'F')) {
                ++pos;
            }
            return emitConst(val);
        }

        if (pos < src.size() && (std::isalpha(static_cast<unsigned char>(src[pos])) || src[pos] == '_')) {
            size_t start = pos;
            while (pos < src.size() && (std::isalnum(static_cast<unsigned char>(src[pos])) || src[pos] == '_' || src[pos] == '.')) {
                ++pos;
            }
            std::string_view ident = src.substr(start, pos - start);
            if (ident.substr(0, 5) == "Math.") {
                ident = ident.substr(5);
            }

            skipWs();
            if (match('(')) {
                std::vector<int> args;
                skipWs();
                if (!match(')')) {
                    for (;;) {
                        args.push_back(parseExpression());
                        skipWs();
                        if (match(',')) continue;
                        if (match(')')) break;
                        throw std::invalid_argument("Custom SDF expression: expected ',' or ')' in function call");
                    }
                }
                return buildCall(ident, args);
            }

            if (ident == "x" || ident == "px") return regX;
            if (ident == "y" || ident == "py") return regY;
            if (ident == "z" || ident == "pz") return regZ;
            if (ident == "pi" || ident == "PI") return emitConst(3.14159265358979323846f);
            if (ident == "tau" || ident == "TAU") return emitConst(6.28318530717958647692f);
            if (ident == "e" || ident == "E") return emitConst(2.71828182845904523536f);

            throw std::invalid_argument("Custom SDF expression: unknown identifier '" + std::string(ident) + "'");
        }

        throw std::invalid_argument("Custom SDF expression: unexpected token at position " + std::to_string(pos));
    }

    int buildCall(std::string_view fn, const std::vector<int>& args) {
        auto requireArgs = [&](size_t n) {
            if (args.size() != n) {
                throw std::invalid_argument("Custom SDF expression: function '" + std::string(fn) +
                                            "' expects " + std::to_string(n) + " arguments");
            }
        };

        if (fn == "sin" || fn == "sinf") {
            requireArgs(1);
            return emit({SdfExprOp::Sin, args[0], -1, -1, 0.0f});
        }
        if (fn == "cos" || fn == "cosf") {
            requireArgs(1);
            return emit({SdfExprOp::Cos, args[0], -1, -1, 0.0f});
        }
        if (fn == "abs" || fn == "fabs" || fn == "fabsf") {
            requireArgs(1);
            return emit({SdfExprOp::Abs, args[0], -1, -1, 0.0f});
        }
        if (fn == "sqrt" || fn == "sqrtf") {
            requireArgs(1);
            return emit({SdfExprOp::Sqrt, args[0], -1, -1, 0.0f});
        }
        if (fn == "floor" || fn == "floorf") {
            requireArgs(1);
            return emit({SdfExprOp::Floor, args[0], -1, -1, 0.0f});
        }
        if (fn == "min" || fn == "fmin" || fn == "fminf") {
            requireArgs(2);
            return emit({SdfExprOp::Min, args[0], args[1], -1, 0.0f});
        }
        if (fn == "max" || fn == "fmax" || fn == "fmaxf") {
            requireArgs(2);
            return emit({SdfExprOp::Max, args[0], args[1], -1, 0.0f});
        }
        if (fn == "clamp") {
            requireArgs(3);
            return emit({SdfExprOp::Clamp, args[0], args[1], args[2], 0.0f});
        }
        if (fn == "fma" || fn == "fmaf") {
            requireArgs(3);
            return emit({SdfExprOp::Fma, args[0], args[1], args[2], 0.0f});
        }
        if (fn == "lerp" || fn == "mix") {
            requireArgs(3);
            int diff = emit({SdfExprOp::Sub, args[1], args[0], -1, 0.0f});
            return emit({SdfExprOp::Fma, diff, args[2], args[0], 0.0f});
        }
        if (fn == "length" || fn == "hypot") {
            if (args.size() == 2) {
                int x2 = emit({SdfExprOp::Mul, args[0], args[0], -1, 0.0f});
                int sum = emit({SdfExprOp::Fma, args[1], args[1], x2, 0.0f});
                return emit({SdfExprOp::Sqrt, sum, -1, -1, 0.0f});
            }
            if (args.size() == 3) {
                int x2 = emit({SdfExprOp::Mul, args[0], args[0], -1, 0.0f});
                int xy2 = emit({SdfExprOp::Fma, args[1], args[1], x2, 0.0f});
                int sum = emit({SdfExprOp::Fma, args[2], args[2], xy2, 0.0f});
                return emit({SdfExprOp::Sqrt, sum, -1, -1, 0.0f});
            }
            throw std::invalid_argument("Custom SDF expression: 'length' expects 2 or 3 arguments");
        }
        if (fn == "pow") {
            requireArgs(2);
            return parseExponentTailFromArgs(args[0], args[1]);
        }
        throw std::invalid_argument("Custom SDF expression: unsupported function '" + std::string(fn) + "'");
    }

    int parseExponentTailFromArgs(int base, int expNode) {
        const auto& expInst = insts[static_cast<size_t>(expNode)];
        if (expInst.op == SdfExprOp::Const) {
            float e = expInst.imm;
            if (e == 2.0f) return emit({SdfExprOp::Mul, base, base, -1, 0.0f});
            if (e == 3.0f) {
                int sq = emit({SdfExprOp::Mul, base, base, -1, 0.0f});
                return emit({SdfExprOp::Mul, sq, base, -1, 0.0f});
            }
            if (e == 4.0f) {
                int sq = emit({SdfExprOp::Mul, base, base, -1, 0.0f});
                return emit({SdfExprOp::Mul, sq, sq, -1, 0.0f});
            }
            if (e == 0.5f) return emit({SdfExprOp::Sqrt, base, -1, -1, 0.0f});
        }
        return emit({SdfExprOp::Mul, base, base, -1, 0.0f});
    }
};

} // namespace

SdfCustomExpr parseCustomSdfExpr(std::string_view expr, float lipschitzBound) {
    ExprParser parser(expr);
    parser.stripJsFunctionHeader();
    int res = parser.parseExpression();
    parser.skipWs();
    if (parser.match(';')) parser.skipWs();
    if (parser.match('}')) parser.skipWs();
    if (parser.pos != parser.src.size()) {
        throw std::invalid_argument("Custom SDF expression: unexpected trailing characters at position " +
                                    std::to_string(parser.pos));
    }

    SdfCustomExpr compiled;
    compiled.source = std::string(expr);
    compiled.instructions = std::move(parser.insts);
    compiled.resultReg = res;
    compiled.lipschitzBound = (std::isfinite(lipschitzBound) && lipschitzBound > 0.0f) ? lipschitzBound : 1.0f;
    return compiled;
}

float evalCustomExprScalar(const SdfCustomExpr& expr, float px, float py, float pz) {
    if (expr.instructions.empty() || expr.resultReg < 0) return 1e9f;
    std::vector<float> regs(expr.instructions.size(), 0.0f);
    for (size_t i = 0; i < expr.instructions.size(); ++i) {
        const SdfExprInst& inst = expr.instructions[i];
        float a = (inst.a >= 0) ? regs[static_cast<size_t>(inst.a)] : 0.0f;
        float b = (inst.b >= 0) ? regs[static_cast<size_t>(inst.b)] : 0.0f;
        float c = (inst.c >= 0) ? regs[static_cast<size_t>(inst.c)] : 0.0f;
        switch (inst.op) {
            case SdfExprOp::VarX:  regs[i] = px; break;
            case SdfExprOp::VarY:  regs[i] = py; break;
            case SdfExprOp::VarZ:  regs[i] = pz; break;
            case SdfExprOp::Const: regs[i] = inst.imm; break;
            case SdfExprOp::Neg:   regs[i] = -a; break;
            case SdfExprOp::Add:   regs[i] = a + b; break;
            case SdfExprOp::Sub:   regs[i] = a - b; break;
            case SdfExprOp::Mul:   regs[i] = a * b; break;
            case SdfExprOp::Div:   regs[i] = (std::fabs(b) > 1e-20f) ? (a / b) : 0.0f; break;
            case SdfExprOp::Sin:   regs[i] = std::sin(a); break;
            case SdfExprOp::Cos:   regs[i] = std::cos(a); break;
            case SdfExprOp::Abs:   regs[i] = std::fabs(a); break;
            case SdfExprOp::Sqrt:  regs[i] = std::sqrt(std::max(a, 0.0f)); break;
            case SdfExprOp::Floor: regs[i] = std::floor(a); break;
            case SdfExprOp::Min:   regs[i] = std::min(a, b); break;
            case SdfExprOp::Max:   regs[i] = std::max(a, b); break;
            case SdfExprOp::Clamp: regs[i] = std::clamp(a, b, c); break;
            case SdfExprOp::Fma:   regs[i] = std::fma(a, b, c); break;
        }
    }
    return regs[static_cast<size_t>(expr.resultReg)];
}

float evalCustomExprGradScalar(const SdfCustomExpr& expr, float px, float py, float pz,
                               float& gx, float& gy, float& gz) {
    if (expr.instructions.empty() || expr.resultReg < 0) {
        gx = gy = gz = 0.0f;
        return 1e9f;
    }
    struct ScalarDual { float v, dx, dy, dz; };
    std::vector<ScalarDual> regs(expr.instructions.size(), {0.0f, 0.0f, 0.0f, 0.0f});
    for (size_t i = 0; i < expr.instructions.size(); ++i) {
        const SdfExprInst& inst = expr.instructions[i];
        ScalarDual a = (inst.a >= 0) ? regs[static_cast<size_t>(inst.a)] : ScalarDual{0, 0, 0, 0};
        ScalarDual b = (inst.b >= 0) ? regs[static_cast<size_t>(inst.b)] : ScalarDual{0, 0, 0, 0};
        ScalarDual c = (inst.c >= 0) ? regs[static_cast<size_t>(inst.c)] : ScalarDual{0, 0, 0, 0};
        switch (inst.op) {
            case SdfExprOp::VarX:  regs[i] = {px, 1.0f, 0.0f, 0.0f}; break;
            case SdfExprOp::VarY:  regs[i] = {py, 0.0f, 1.0f, 0.0f}; break;
            case SdfExprOp::VarZ:  regs[i] = {pz, 0.0f, 0.0f, 1.0f}; break;
            case SdfExprOp::Const: regs[i] = {inst.imm, 0.0f, 0.0f, 0.0f}; break;
            case SdfExprOp::Neg:   regs[i] = {-a.v, -a.dx, -a.dy, -a.dz}; break;
            case SdfExprOp::Add:   regs[i] = {a.v + b.v, a.dx + b.dx, a.dy + b.dy, a.dz + b.dz}; break;
            case SdfExprOp::Sub:   regs[i] = {a.v - b.v, a.dx - b.dx, a.dy - b.dy, a.dz - b.dz}; break;
            case SdfExprOp::Mul:
                regs[i] = {a.v * b.v, a.v * b.dx + b.v * a.dx, a.v * b.dy + b.v * a.dy, a.v * b.dz + b.v * a.dz};
                break;
            case SdfExprOp::Div: {
                float inv = (std::fabs(b.v) > 1e-20f) ? (1.0f / b.v) : 0.0f;
                float val = a.v * inv;
                regs[i] = {val, (a.dx - val * b.dx) * inv, (a.dy - val * b.dy) * inv, (a.dz - val * b.dz) * inv};
                break;
            }
            case SdfExprOp::Sin: {
                float s = std::sin(a.v), co = std::cos(a.v);
                regs[i] = {s, co * a.dx, co * a.dy, co * a.dz};
                break;
            }
            case SdfExprOp::Cos: {
                float s = std::sin(a.v), co = std::cos(a.v);
                regs[i] = {co, -s * a.dx, -s * a.dy, -s * a.dz};
                break;
            }
            case SdfExprOp::Abs: {
                float s = (a.v >= 0.0f) ? 1.0f : -1.0f;
                regs[i] = {std::fabs(a.v), s * a.dx, s * a.dy, s * a.dz};
                break;
            }
            case SdfExprOp::Sqrt: {
                float r = std::sqrt(std::max(a.v, 0.0f));
                float inv2r = (r > 1e-12f) ? (0.5f / r) : 0.0f;
                regs[i] = {r, inv2r * a.dx, inv2r * a.dy, inv2r * a.dz};
                break;
            }
            case SdfExprOp::Floor:
                regs[i] = {std::floor(a.v), 0.0f, 0.0f, 0.0f};
                break;
            case SdfExprOp::Min:
                regs[i] = (a.v <= b.v) ? a : b;
                break;
            case SdfExprOp::Max:
                regs[i] = (a.v >= b.v) ? a : b;
                break;
            case SdfExprOp::Clamp: {
                ScalarDual mn = (a.v <= c.v) ? a : c;
                regs[i] = (mn.v >= b.v) ? mn : b;
                break;
            }
            case SdfExprOp::Fma:
                regs[i] = {std::fma(a.v, b.v, c.v),
                           a.v * b.dx + b.v * a.dx + c.dx,
                           a.v * b.dy + b.v * a.dy + c.dy,
                           a.v * b.dz + b.v * a.dz + c.dz};
                break;
        }
    }
    const ScalarDual& out = regs[static_cast<size_t>(expr.resultReg)];
    gx = out.dx;
    gy = out.dy;
    gz = out.dz;
    return out.v;
}

#if BROMESH_HAS_BRASS_JIT
using namespace brass;
using namespace brass::codegen;

// Branchless Cody-Waite range reduction + degree-9 minimax polynomial for sin(x) and cos(x).
// On GPU targets (SPIR-V / PTX), emits native GLSL.std.450 Sin/Cos / PTX sin.approx/cos.approx.
Value* emitSinF32(SdfKernelBuilder& kb, Value* x, bool gpuIntrinsic) {
    if (gpuIntrinsic) {
        if (kb.builder().current_module() && !kb.builder().current_module()->has_external_symbol("sinf")) {
            kb.builder().current_module()->add_external_symbol("sinf");
        }
        return kb.builder().build_call("sinf", Type::f32(), {x});
    }
    Value* invPi = kb.const_f32(0.31830988618379067f);
    Value* k = kb.builder().build_round_f32(kb.mul(x, invPi));
    Value* negPiHi = kb.const_f32(-3.140625f);
    Value* negPiLo = kb.const_f32(-0.00096765358979323846f);
    Value* r = kb.fma(k, negPiLo, kb.fma(k, negPiHi, x));
    Value* r2 = kb.mul(r, r);
    Value* r3 = kb.mul(r2, r);

    Value* c9 = kb.const_f32(2.7525562e-6f);
    Value* c7 = kb.const_f32(-0.00019840874f);
    Value* c5 = kb.const_f32(0.0083333310f);
    Value* c3 = kb.const_f32(-0.16666667f);

    Value* p = kb.fma(r2, c9, c7);
    p = kb.fma(r2, p, c5);
    p = kb.fma(r2, p, c3);
    p = kb.fma(r3, p, r);

    Value* ki = kb.builder().build_fptosi_i32_f32(k);
    Value* signBit = kb.builder().build_shl(kb.builder().build_and(ki, kb.const_i32(1)), kb.const_i32(31));
    Value* pBits = kb.builder().build_bitcast_i32_f32(p);
    return kb.builder().build_bitcast_f32_i32(kb.builder().build_xor(pBits, signBit));
}

Value* emitCosF32(SdfKernelBuilder& kb, Value* x, bool gpuIntrinsic) {
    if (gpuIntrinsic) {
        if (kb.builder().current_module() && !kb.builder().current_module()->has_external_symbol("cosf")) {
            kb.builder().current_module()->add_external_symbol("cosf");
        }
        return kb.builder().build_call("cosf", Type::f32(), {x});
    }
    return emitSinF32(kb, kb.add(x, kb.const_f32(1.5707963267948966f)), false);
}

Value* lowerCustomExprToMir(SdfKernelBuilder& kb, const SdfCustomExpr& expr,
                            Value* px, Value* py, Value* pz, bool gpuTarget) {
    if (expr.instructions.empty() || expr.resultReg < 0) {
        return kb.const_f32(1e9f);
    }
    std::vector<Value*> regs(expr.instructions.size(), nullptr);
    for (size_t i = 0; i < expr.instructions.size(); ++i) {
        const SdfExprInst& inst = expr.instructions[i];
        Value* a = (inst.a >= 0) ? regs[static_cast<size_t>(inst.a)] : nullptr;
        Value* b = (inst.b >= 0) ? regs[static_cast<size_t>(inst.b)] : nullptr;
        Value* c = (inst.c >= 0) ? regs[static_cast<size_t>(inst.c)] : nullptr;
        switch (inst.op) {
            case SdfExprOp::VarX:  regs[i] = px; break;
            case SdfExprOp::VarY:  regs[i] = py; break;
            case SdfExprOp::VarZ:  regs[i] = pz; break;
            case SdfExprOp::Const: regs[i] = kb.const_f32(inst.imm); break;
            case SdfExprOp::Neg:   regs[i] = kb.builder().build_neg(a); break;
            case SdfExprOp::Add:   regs[i] = kb.add(a, b); break;
            case SdfExprOp::Sub:   regs[i] = kb.sub(a, b); break;
            case SdfExprOp::Mul:   regs[i] = kb.mul(a, b); break;
            case SdfExprOp::Div:   regs[i] = kb.div(a, b); break;
            case SdfExprOp::Sin:   regs[i] = emitSinF32(kb, a, gpuTarget); break;
            case SdfExprOp::Cos:   regs[i] = emitCosF32(kb, a, gpuTarget); break;
            case SdfExprOp::Abs:   regs[i] = kb.abs_f32(a); break;
            case SdfExprOp::Sqrt:  regs[i] = kb.sqrt_f32(kb.max_f32(a, kb.const_f32(0.0f))); break;
            case SdfExprOp::Floor: regs[i] = kb.floor_f32(a); break;
            case SdfExprOp::Min:   regs[i] = kb.min_f32(a, b); break;
            case SdfExprOp::Max:   regs[i] = kb.max_f32(a, b); break;
            case SdfExprOp::Clamp: regs[i] = kb.clamp_f32(a, b, c); break;
            case SdfExprOp::Fma:   regs[i] = kb.fma(a, b, c); break;
        }
    }
    return regs[static_cast<size_t>(expr.resultReg)];
}
#endif

} // namespace bromesh
