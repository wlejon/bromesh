#pragma once

#include "bromesh/isosurface/jit/sdf_compiler.h"

namespace bromesh {

#if BROMESH_HAS_BRASS_JIT
struct DualF32 {
    brass::Value* val = nullptr;
    brass::Value* dx = nullptr;
    brass::Value* dy = nullptr;
    brass::Value* dz = nullptr;
};

struct MirInterval {
    brass::Value* lo = nullptr;
    brass::Value* hi = nullptr;
};

brass::Value* emitSinF32(brass::codegen::SdfKernelBuilder& kb, brass::Value* x, bool gpuIntrinsic = false);
brass::Value* emitCosF32(brass::codegen::SdfKernelBuilder& kb, brass::Value* x, bool gpuIntrinsic = false);

brass::Value* lowerNoise3d(brass::codegen::SdfKernelBuilder& kb, brass::Value* px, brass::Value* py, brass::Value* pz);

brass::Value* lowerToMir(brass::codegen::SdfKernelBuilder& kb, const SdfGraph& graph, int nodeId,
                         brass::Value* px, brass::Value* py, brass::Value* pz, bool gpuTarget = false);

brass::Value* lowerCustomExprToMir(brass::codegen::SdfKernelBuilder& kb, const SdfCustomExpr& expr,
                                   brass::Value* px, brass::Value* py, brass::Value* pz, bool gpuTarget = false);

DualF32 lowerDualToMir(brass::codegen::SdfKernelBuilder& kb, const SdfGraph& graph, int nodeId,
                       DualF32 px, DualF32 py, DualF32 pz);

DualF32 lowerCustomExprDualToMir(brass::codegen::SdfKernelBuilder& kb, const SdfCustomExpr& expr,
                                 DualF32 px, DualF32 py, DualF32 pz);

brass::Function* buildSdfGradFunction(brass::Module& mod, const SdfGraph& graph, int rootNode);

MirInterval lowerIntervalToMir(brass::codegen::SdfKernelBuilder& kb, const SdfGraph& graph, int nodeId,
                               MirInterval bx, MirInterval by, MirInterval bz);

MirInterval lowerCustomExprIntervalToMir(brass::codegen::SdfKernelBuilder& kb, const SdfCustomExpr& expr,
                                         MirInterval bx, MirInterval by, MirInterval bz);

brass::Function* buildSdfBlockBoundsFunction(brass::Module& mod, const SdfGraph& graph, int rootNode);
#endif

float evalScalarNode(const SdfGraph& graph, int nodeId, float px, float py, float pz);
float evalCustomExprScalar(const SdfCustomExpr& expr, float px, float py, float pz);

float evalGradScalarNode(const SdfGraph& graph, int nodeId, float px, float py, float pz,
                         float& gx, float& gy, float& gz);
float evalCustomExprGradScalar(const SdfCustomExpr& expr, float px, float py, float pz,
                               float& gx, float& gy, float& gz);

SdfInterval evalCustomExprInterval(const SdfCustomExpr& expr, SdfInterval bx, SdfInterval by, SdfInterval bz);

} // namespace bromesh
