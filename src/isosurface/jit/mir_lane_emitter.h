#pragma once

#include "bromesh/isosurface/jit/sdf_node.h"

#if BROMESH_HAS_BRASS_JIT
#include <brass/codegen/sdf_builder.hpp>

#include <cstdint>

namespace bromesh::mir_detail {

/// Single-source Brass MIR emitter parameterized over SIMD vector lane width.
/// When lanes == 1 (GPU SIMT or CPU scalar tail), emits scalar f32/i32 MIR.
/// When lanes == 4 or 8 (CPU SIMD), emits f32x4/f32x8 and i32x4/i32x8 MIR.
class MirLaneEmitter {
public:
    MirLaneEmitter(brass::codegen::SdfKernelBuilder& kb, uint32_t lanes, bool gpuMode)
        : kb_(kb), b_(kb.builder()), lanes_(lanes), gpuMode_(gpuMode) {}

    uint32_t lanes() const noexcept { return lanes_; }
    bool isGpu() const noexcept { return gpuMode_; }
    brass::codegen::SdfKernelBuilder& kb() noexcept { return kb_; }
    brass::Builder& b() noexcept { return b_; }

    brass::Type f32Type() const noexcept;
    brass::Type i32Type() const noexcept;

    brass::Value* broadcastF32(brass::Value* scalar);
    brass::Value* cF32(float v);
    brass::Value* cI32(int32_t v);

    brass::Value* add(brass::Value* a, brass::Value* b);
    brass::Value* sub(brass::Value* a, brass::Value* b);
    brass::Value* mul(brass::Value* a, brass::Value* b);
    brass::Value* fma(brass::Value* a, brass::Value* b, brass::Value* c);
    brass::Value* neg(brass::Value* a);
    brass::Value* fmin(brass::Value* a, brass::Value* b);
    brass::Value* fmax(brass::Value* a, brass::Value* b);
    brass::Value* sqrt(brass::Value* a);
    brass::Value* floor(brass::Value* a);
    brass::Value* round(brass::Value* a);
    brass::Value* abs(brass::Value* a);

    brass::Value* rawDiv(brass::Value* a, brass::Value* b);
    brass::Value* safeDiv(brass::Value* a, brass::Value* b);

    brass::Value* clamp(brass::Value* x, brass::Value* lo, brass::Value* hi);
    brass::Value* lerp(brass::Value* a, brass::Value* b, brass::Value* t);
    brass::Value* length2(brass::Value* x, brass::Value* y);
    brass::Value* length3(brass::Value* x, brass::Value* y, brass::Value* z);
    brass::Value* selectGt(brass::Value* lhs, brass::Value* rhs, brass::Value* tVal, brass::Value* fVal);

    brass::Value* fptosiI32(brass::Value* a);
    brass::Value* sitofpF32(brass::Value* a);
    brass::Value* bitcastI32(brass::Value* a);
    brass::Value* bitcastF32(brass::Value* a);

    brass::Value* ixor(brass::Value* a, brass::Value* b);
    brass::Value* iand(brass::Value* a, brass::Value* b);
    brass::Value* ishl(brass::Value* a, brass::Value* b);
    brass::Value* ilshr(brass::Value* a, brass::Value* b);

    brass::Value* sinF32(brass::Value* x);
    brass::Value* cosF32(brass::Value* x);

    brass::Value* loadLaneF32(brass::Value* basePtr, brass::Value* elemIdxI32);
    void storeLaneF32(brass::Value* basePtr, brass::Value* elemIdxI32, brass::Value* val);
    brass::Value* loadScalarBroadcastF32(brass::Value* ptr, int32_t byteOffset);

private:
    brass::codegen::SdfKernelBuilder& kb_;
    brass::Builder& b_;
    uint32_t lanes_;
    bool gpuMode_;
};

brass::Value* buildLaneRampF32(MirLaneEmitter& em);

brass::Value* emitSdfGraphSingleSource(MirLaneEmitter& em, const SdfGraph& graph, int nodeId,
                                       brass::Value* px, brass::Value* py, brass::Value* pz);

uint32_t detectHostSimdLanes();

brass::Function* buildSdfGridSimdFunction(brass::Module& mod, const SdfGraph& graph, int rootNode,
                                          uint32_t lanes, const char* fnName = "eval_sdf_grid_simd");

} // namespace bromesh::mir_detail
#endif
