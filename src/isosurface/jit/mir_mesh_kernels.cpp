#include "bromesh/isosurface/jit/mir_mesh_kernels.h"
#include "mir_lane_emitter.h"

#if BROMESH_HAS_BRASS_JIT
#include <brass/core/diagnostics.hpp>
#include <brass/gpu/vulkan_driver.hpp>
#include <brass/mir/gpu_opt.hpp>
#include <brass/mir/gpu_portability.hpp>
#endif

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <iostream>

namespace bromesh {

namespace {

using GridSimdCpuFn = void (*)(float* outField, int32_t dimX, int32_t dimY, int32_t dimZ,
                               float minX, float minY, float minZ,
                               float stepX, float stepY, float stepZ);

using VertTransformCpuFn = void (*)(float* outSoa, const float* inSoa, const float* matrices25,
                                    int32_t numVerts, int32_t strideElems);

constexpr uint32_t kGpuWorkgroupSize = 64;

void computeGridOriginAndStep(int dimX, int dimY, int dimZ, const bromath::AABB3& bounds,
                              float outOrigin[3], float outStep[3]) {
    outOrigin[0] = bounds.min.x;
    outOrigin[1] = bounds.min.y;
    outOrigin[2] = bounds.min.z;
    outStep[0] = (dimX > 1) ? (bounds.max.x - bounds.min.x) / static_cast<float>(dimX - 1) : 1.0f;
    outStep[1] = (dimY > 1) ? (bounds.max.y - bounds.min.y) / static_cast<float>(dimY - 1) : 1.0f;
    outStep[2] = (dimZ > 1) ? (bounds.max.z - bounds.min.z) / static_cast<float>(dimZ - 1) : 1.0f;
}

int32_t alignedSoaStride(size_t numVerts) {
    // Pad channel stride to 8 floats (32 bytes) so AVX2 loads stay within channel bounds.
    const size_t padded = ((numVerts + 7u) / 8u) * 8u;
    return static_cast<int32_t>(std::max<size_t>(8u, padded));
}

void packMeshToSoa(const MeshData& mesh, std::vector<float>& outSoa, int32_t strideElems) {
    const size_t numVerts = mesh.positions.size() / 3;
    const bool hasNormals = mesh.normals.size() == mesh.positions.size();
    outSoa.assign(static_cast<size_t>(strideElems) * 6u, 0.0f);
    float* px = outSoa.data() + 0 * strideElems;
    float* py = outSoa.data() + 1 * strideElems;
    float* pz = outSoa.data() + 2 * strideElems;
    float* nx = outSoa.data() + 3 * strideElems;
    float* ny = outSoa.data() + 4 * strideElems;
    float* nz = outSoa.data() + 5 * strideElems;

    for (size_t i = 0; i < numVerts; ++i) {
        px[i] = mesh.positions[3 * i + 0];
        py[i] = mesh.positions[3 * i + 1];
        pz[i] = mesh.positions[3 * i + 2];
        nx[i] = hasNormals ? mesh.normals[3 * i + 0] : 0.0f;
        ny[i] = hasNormals ? mesh.normals[3 * i + 1] : 1.0f;
        nz[i] = hasNormals ? mesh.normals[3 * i + 2] : 0.0f;
    }
}

void unpackSoaToMesh(const std::vector<float>& soa, int32_t strideElems, size_t numVerts, MeshData& dst) {
    dst.positions.resize(numVerts * 3);
    dst.normals.resize(numVerts * 3);
    const float* px = soa.data() + 0 * strideElems;
    const float* py = soa.data() + 1 * strideElems;
    const float* pz = soa.data() + 2 * strideElems;
    const float* nx = soa.data() + 3 * strideElems;
    const float* ny = soa.data() + 4 * strideElems;
    const float* nz = soa.data() + 5 * strideElems;

    for (size_t i = 0; i < numVerts; ++i) {
        dst.positions[3 * i + 0] = px[i];
        dst.positions[3 * i + 1] = py[i];
        dst.positions[3 * i + 2] = pz[i];
        dst.normals[3 * i + 0] = nx[i];
        dst.normals[3 * i + 1] = ny[i];
        dst.normals[3 * i + 2] = nz[i];
    }
}

template <typename ExecFn>
bool transformMeshViaSoa(const MeshData& src, MeshData& dst, const float mat4[16], ExecFn&& exec) {
    dst = src;
    const size_t numVerts = src.positions.size() / 3;
    if (numVerts == 0) {
        return true;
    }
    const int32_t stride = alignedSoaStride(numVerts);
    std::vector<float> inSoa;
    std::vector<float> outSoa(static_cast<size_t>(stride) * 6u, 0.0f);
    packMeshToSoa(src, inSoa, stride);

    float mats[25]{};
    std::memcpy(mats, mat4, 16u * sizeof(float));
    computeNormalMatrix3x3(mat4, mats + 16);

    if (!exec(inSoa.data(), outSoa.data(), mats, static_cast<int32_t>(numVerts), stride)) {
        return false;
    }
    unpackSoaToMesh(outSoa, stride, numVerts, dst);
    return true;
}

} // namespace

void computeNormalMatrix3x3(const float mat4[16], float outNormMat3[9]) {
    const float a00 = mat4[0], a10 = mat4[1], a20 = mat4[2];
    const float a01 = mat4[4], a11 = mat4[5], a21 = mat4[6];
    const float a02 = mat4[8], a12 = mat4[9], a22 = mat4[10];

    const float c00 = a11 * a22 - a12 * a21;
    const float c01 = -(a10 * a22 - a12 * a20);
    const float c02 = a10 * a21 - a11 * a20;

    const float c10 = -(a01 * a22 - a02 * a21);
    const float c11 = a00 * a22 - a02 * a20;
    const float c12 = -(a00 * a21 - a01 * a20);

    const float c20 = a01 * a12 - a02 * a11;
    const float c21 = -(a00 * a12 - a02 * a10);
    const float c22 = a00 * a11 - a01 * a10;

    const float det = a00 * c00 + a01 * c01 + a02 * c02;
    if (std::fabs(det) < 1e-12f) {
        outNormMat3[0] = 1.0f; outNormMat3[1] = 0.0f; outNormMat3[2] = 0.0f;
        outNormMat3[3] = 0.0f; outNormMat3[4] = 1.0f; outNormMat3[5] = 0.0f;
        outNormMat3[6] = 0.0f; outNormMat3[7] = 0.0f; outNormMat3[8] = 1.0f;
        return;
    }
    const float invDet = 1.0f / det;
    outNormMat3[0] = c00 * invDet; outNormMat3[1] = c10 * invDet; outNormMat3[2] = c20 * invDet;
    outNormMat3[3] = c01 * invDet; outNormMat3[4] = c11 * invDet; outNormMat3[5] = c21 * invDet;
    outNormMat3[6] = c02 * invDet; outNormMat3[7] = c12 * invDet; outNormMat3[8] = c22 * invDet;
}

void transformMeshReference(const MeshData& src, MeshData& dst, const float mat4[16]) {
    dst = src;
    const size_t numVerts = src.positions.size() / 3;
    const bool hasNormals = src.normals.size() == src.positions.size();
    dst.normals.resize(src.positions.size());

    float nmat[9]{};
    computeNormalMatrix3x3(mat4, nmat);

    for (size_t i = 0; i < numVerts; ++i) {
        const float px = src.positions[3 * i + 0];
        const float py = src.positions[3 * i + 1];
        const float pz = src.positions[3 * i + 2];

        dst.positions[3 * i + 0] = mat4[0] * px + mat4[4] * py + mat4[8] * pz + mat4[12];
        dst.positions[3 * i + 1] = mat4[1] * px + mat4[5] * py + mat4[9] * pz + mat4[13];
        dst.positions[3 * i + 2] = mat4[2] * px + mat4[6] * py + mat4[10] * pz + mat4[14];

        const float nx = hasNormals ? src.normals[3 * i + 0] : 0.0f;
        const float ny = hasNormals ? src.normals[3 * i + 1] : 1.0f;
        const float nz = hasNormals ? src.normals[3 * i + 2] : 0.0f;

        const float tx = nmat[0] * nx + nmat[3] * ny + nmat[6] * nz;
        const float ty = nmat[1] * nx + nmat[4] * ny + nmat[7] * nz;
        const float tz = nmat[2] * nx + nmat[5] * ny + nmat[8] * nz;
        const float len = std::sqrt(tx * tx + ty * ty + tz * tz);
        if (len > 1e-7f) {
            const float inv = 1.0f / len;
            dst.normals[3 * i + 0] = tx * inv;
            dst.normals[3 * i + 1] = ty * inv;
            dst.normals[3 * i + 2] = tz * inv;
        } else {
            dst.normals[3 * i + 0] = 0.0f;
            dst.normals[3 * i + 1] = 1.0f;
            dst.normals[3 * i + 2] = 0.0f;
        }
    }
}

#if BROMESH_HAS_BRASS_JIT
namespace {

using namespace brass;
using namespace brass::codegen;
using mir_detail::MirLaneEmitter;
using mir_detail::buildLaneRampF32;
using mir_detail::emitSdfGraphSingleSource;

uint32_t detectHostSimdLanesImpl() {
    if (Target::host().is_x64()) {
        return Target::host_has_avx() ? 8u : 4u;
    }
    if (Target::host().is_aarch64()) {
        return 4u;
    }
    return 1u;
}

KernelOptions makeSimdJitOptions(uint32_t lanes) {
    KernelOptions opts;
    opts.opt_level = 2;
    opts.enable_optimizations = true;
    opts.enable_sroa = true;
    opts.enable_bce = true;
    opts.enable_smt_guard_elim = true;
    opts.enable_superopt = false;
    opts.enable_avx2 = (lanes == 8u);
    opts.enable_fma = true;
    opts.enable_vectorize = false;
    opts.enable_slp = false;
    opts.enable_unroll = false;
    opts.enable_parallel = false;
    return opts;
}

Function* buildGridSdfCpuSimdFunction(Module& mod, const SdfGraph& graph, int root, uint32_t lanes) {
    return mir_detail::buildSdfGridSimdFunction(mod, graph, root, lanes, "eval_sdf_grid_simd");
}

Function* buildGridSdfGpuFunction(Module& mod, const SdfGraph& graph, int root) {
    Function* fn = mod.create_function("eval_sdf_grid_gpu", Type::void_type(), {
        Type::ptr(), Type::i32(), Type::i32(), Type::i32(),
        Type::f32(), Type::f32(), Type::f32(),
        Type::f32(), Type::f32(), Type::f32()
    });
    SdfKernelBuilder kb(mod, fn);
    Builder& b = kb.builder();
    BasicBlock* entry = b.append_block("entry");
    Value* outField = b.add_block_param(entry, Type::ptr());
    Value* dimX = b.add_block_param(entry, Type::i32());
    Value* dimY = b.add_block_param(entry, Type::i32());
    Value* dimZ = b.add_block_param(entry, Type::i32());
    Value* minX = b.add_block_param(entry, Type::f32());
    Value* minY = b.add_block_param(entry, Type::f32());
    Value* minZ = b.add_block_param(entry, Type::f32());
    Value* stepX = b.add_block_param(entry, Type::f32());
    Value* stepY = b.add_block_param(entry, Type::f32());
    Value* stepZ = b.add_block_param(entry, Type::f32());

    kb.position_at_end(entry);
    Value* idx = b.build_call("gpu_global_tid_x", Type::i32());
    Value* strideZ = b.build_mul(dimX, dimY);
    Value* total = b.build_mul(strideZ, dimZ);
    Value* inBounds = b.build_ult(idx, total);

    kb.if_then(inBounds, [&]() {
        Value* iz = b.build_sdiv(idx, strideZ);
        Value* rem = b.build_sub(idx, b.build_mul(iz, strideZ));
        Value* iy = b.build_sdiv(rem, dimX);
        Value* ix = b.build_sub(rem, b.build_mul(iy, dimX));

        Value* px = b.build_fma(kb.to_f32(ix), stepX, minX);
        Value* py = b.build_fma(kb.to_f32(iy), stepY, minY);
        Value* pz = b.build_fma(kb.to_f32(iz), stepZ, minZ);

        MirLaneEmitter emGpu(kb, 1, true);
        Value* sdfVal = emitSdfGraphSingleSource(emGpu, graph, root, px, py, pz);
        emGpu.storeLaneF32(outField, idx, sdfVal);
    });

    b.build_ret_void();
    fn->rebuild_cfg_predecessors();
    return fn;
}

void emitVertexTransformAndNormalizeSingleSource(MirLaneEmitter& em,
                                                 Value* outSoa, Value* inSoa,
                                                 Value* mats, Value* strideElemsI32,
                                                 Value* elemIdxI32) {
    Value* idx1 = em.b().build_add(elemIdxI32, strideElemsI32);
    Value* twoStride = em.b().build_mul(strideElemsI32, em.kb().const_i32(2));
    Value* idx2 = em.b().build_add(elemIdxI32, twoStride);
    Value* threeStride = em.b().build_mul(strideElemsI32, em.kb().const_i32(3));
    Value* idx3 = em.b().build_add(elemIdxI32, threeStride);
    Value* idx4 = em.b().build_add(idx1, threeStride);
    Value* idx5 = em.b().build_add(idx2, threeStride);

    Value* px = em.loadLaneF32(inSoa, elemIdxI32);
    Value* py = em.loadLaneF32(inSoa, idx1);
    Value* pz = em.loadLaneF32(inSoa, idx2);

    auto mLoad = [&](int idx) { return em.loadScalarBroadcastF32(mats, idx * 4); };
    Value* ox = em.fma(mLoad(0), px, em.fma(mLoad(4), py, em.fma(mLoad(8), pz, mLoad(12))));
    Value* oy = em.fma(mLoad(1), px, em.fma(mLoad(5), py, em.fma(mLoad(9), pz, mLoad(13))));
    Value* oz = em.fma(mLoad(2), px, em.fma(mLoad(6), py, em.fma(mLoad(10), pz, mLoad(14))));

    em.storeLaneF32(outSoa, elemIdxI32, ox);
    em.storeLaneF32(outSoa, idx1, oy);
    em.storeLaneF32(outSoa, idx2, oz);

    Value* nx = em.loadLaneF32(inSoa, idx3);
    Value* ny = em.loadLaneF32(inSoa, idx4);
    Value* nz = em.loadLaneF32(inSoa, idx5);

    auto nLoad = [&](int idx) { return em.loadScalarBroadcastF32(mats, (16 + idx) * 4); };
    Value* tx = em.fma(nLoad(0), nx, em.fma(nLoad(3), ny, em.mul(nLoad(6), nz)));
    Value* ty = em.fma(nLoad(1), nx, em.fma(nLoad(4), ny, em.mul(nLoad(7), nz)));
    Value* tz = em.fma(nLoad(2), nx, em.fma(nLoad(5), ny, em.mul(nLoad(8), nz)));

    Value* len = em.length3(tx, ty, tz);
    Value* eps = em.cF32(1e-7f);
    Value* zero = em.cF32(0.0f);
    Value* one = em.cF32(1.0f);
    Value* safeLen = em.selectGt(len, eps, len, one);
    Value* invLen = em.rawDiv(one, safeLen);

    Value* onx = em.selectGt(len, eps, em.mul(tx, invLen), zero);
    Value* ony = em.selectGt(len, eps, em.mul(ty, invLen), one);
    Value* onz = em.selectGt(len, eps, em.mul(tz, invLen), zero);

    em.storeLaneF32(outSoa, idx3, onx);
    em.storeLaneF32(outSoa, idx4, ony);
    em.storeLaneF32(outSoa, idx5, onz);
}

Function* buildVertexTransformCpuSimdFunction(Module& mod, uint32_t lanes) {
    Function* fn = mod.create_function("transform_mesh_simd", Type::void_type(), {
        Type::ptr(), Type::ptr(), Type::ptr(), Type::i32(), Type::i32()
    });
    SdfKernelBuilder kb(mod, fn);
    Builder& b = kb.builder();
    BasicBlock* entry = b.append_block("entry");
    Value* outSoa = b.add_block_param(entry, Type::ptr());
    Value* inSoa = b.add_block_param(entry, Type::ptr());
    Value* mats = b.add_block_param(entry, Type::ptr());
    Value* numVerts = b.add_block_param(entry, Type::i32());
    Value* strideElems = b.add_block_param(entry, Type::i32());

    kb.position_at_end(entry);
    MirLaneEmitter emSimd(kb, lanes, false);
    MirLaneEmitter emScalar(kb, 1, false);

    Value* zero = kb.const_i32(0);
    Value* one = kb.const_i32(1);
    Value* lanesI32 = kb.const_i32(static_cast<int32_t>(lanes));
    Value* vecLimit = (lanes > 1)
                          ? b.build_sub(numVerts, kb.const_i32(static_cast<int32_t>(lanes - 1)))
                          : zero;

    if (lanes > 1) {
        kb.for_range(zero, vecLimit, lanesI32, [&](Value* idx) {
            emitVertexTransformAndNormalizeSingleSource(emSimd, outSoa, inSoa, mats, strideElems, idx);
        });
    }

    Value* remStart = (lanes > 1) ? b.build_mul(b.build_sdiv(numVerts, lanesI32), lanesI32) : zero;
    kb.for_range(remStart, numVerts, one, [&](Value* idx) {
        emitVertexTransformAndNormalizeSingleSource(emScalar, outSoa, inSoa, mats, strideElems, idx);
    });

    b.build_ret_void();
    fn->rebuild_cfg_predecessors();
    return fn;
}

Function* buildVertexTransformGpuFunction(Module& mod) {
    Function* fn = mod.create_function("transform_mesh_gpu", Type::void_type(), {
        Type::ptr(), Type::ptr(), Type::ptr(), Type::i32(), Type::i32()
    });
    SdfKernelBuilder kb(mod, fn);
    Builder& b = kb.builder();
    BasicBlock* entry = b.append_block("entry");
    Value* outSoa = b.add_block_param(entry, Type::ptr());
    Value* inSoa = b.add_block_param(entry, Type::ptr());
    Value* mats = b.add_block_param(entry, Type::ptr());
    Value* numVerts = b.add_block_param(entry, Type::i32());
    Value* strideElems = b.add_block_param(entry, Type::i32());

    kb.position_at_end(entry);
    Value* idx = b.build_call("gpu_global_tid_x", Type::i32());
    Value* inBounds = b.build_ult(idx, numVerts);

    kb.if_then(inBounds, [&]() {
        MirLaneEmitter emGpu(kb, 1, true);
        emitVertexTransformAndNormalizeSingleSource(emGpu, outSoa, inSoa, mats, strideElems, idx);
    });

    b.build_ret_void();
    fn->rebuild_cfg_predecessors();
    return fn;
}

void lowerPortableGpuTargets(Module& mod, Function& fn, MirGpuArtifacts& out,
                             brass::target::SpirvKernel& outSpvKernel) {
    GpuOptOptions gpuOpts;
    gpuOpts.opt_level = 2;
    gpuOpts.enable_fma = false;
    optimize_gpu_function(fn, mod, gpuOpts);

    DiagnosticReporter diag;
    GpuPortabilityOptions portOpts;
    out.portabilityOk = verify_gpu_portability(fn, portOpts, &diag);
    for (const auto& d : diag.diagnostics()) {
        out.portabilityDiagnostics += d.message + "\n";
    }
    if (!out.portabilityOk) {
        return;
    }

    target::SpirvOptions spvOpts;
    spvOpts.stage = target::SpirvStage::Compute;
    spvOpts.local_size_x = kGpuWorkgroupSize;
    outSpvKernel = target::SpirvTarget::compile(fn, spvOpts);
    out.spirvWords = outSpvKernel.words;
    out.spirvEntry = outSpvKernel.entry;

    out.ptxSource = target::PtxTarget::emit_function(fn);

    target::MslOptions mslOpts;
    mslOpts.workgroup_x = kGpuWorkgroupSize;
    target::MslKernel mslKernel = target::MslTarget::compile(fn, mslOpts);
    out.mslSource = std::move(mslKernel.source);
    out.mslEntry = std::move(mslKernel.entry_name);
    out.workgroupX = mslKernel.workgroup_x;
}

} // namespace

namespace mir_detail {

uint32_t detectHostSimdLanes() {
    return detectHostSimdLanesImpl();
}

Function* buildSdfGridSimdFunction(Module& mod, const SdfGraph& graph, int rootNode,
                                   uint32_t lanes, const char* fnName) {
    const int root = (rootNode >= 0) ? rootNode : graph.root();
    Function* fn = mod.create_function(fnName ? fnName : "eval_sdf_grid_simd", Type::void_type(), {
        Type::ptr(), Type::i32(), Type::i32(), Type::i32(),
        Type::f32(), Type::f32(), Type::f32(),
        Type::f32(), Type::f32(), Type::f32()
    });
    SdfKernelBuilder kb(mod, fn);
    Builder& b = kb.builder();
    BasicBlock* entry = b.append_block("entry");
    Value* outField = b.add_block_param(entry, Type::ptr());
    Value* dimX = b.add_block_param(entry, Type::i32());
    Value* dimY = b.add_block_param(entry, Type::i32());
    Value* dimZ = b.add_block_param(entry, Type::i32());
    Value* minX = b.add_block_param(entry, Type::f32());
    Value* minY = b.add_block_param(entry, Type::f32());
    Value* minZ = b.add_block_param(entry, Type::f32());
    Value* stepX = b.add_block_param(entry, Type::f32());
    Value* stepY = b.add_block_param(entry, Type::f32());
    Value* stepZ = b.add_block_param(entry, Type::f32());

    kb.position_at_end(entry);
    MirLaneEmitter emSimd(kb, lanes, false);
    MirLaneEmitter emScalar(kb, 1, false);

    Value* zero = kb.const_i32(0);
    Value* one = kb.const_i32(1);
    Value* lanesI32 = kb.const_i32(static_cast<int32_t>(lanes));
    Value* vecLimit = (lanes > 1)
                          ? b.build_sub(dimX, kb.const_i32(static_cast<int32_t>(lanes - 1)))
                          : zero;
    Value* laneRamp = (lanes > 1) ? buildLaneRampF32(emSimd) : nullptr;
    Value* vMinX = (lanes > 1) ? emSimd.broadcastF32(minX) : nullptr;
    Value* vStepX = (lanes > 1) ? emSimd.broadcastF32(stepX) : nullptr;

    kb.for_range(zero, dimZ, one, [&](Value* iz) {
        Value* pzScalar = b.build_fma(kb.to_f32(iz), stepZ, minZ);
        Value* pzVec = (lanes > 1) ? emSimd.broadcastF32(pzScalar) : nullptr;
        Value* rowZ = b.build_mul(iz, dimY);

        kb.for_range(zero, dimY, one, [&](Value* iy) {
            Value* pyScalar = b.build_fma(kb.to_f32(iy), stepY, minY);
            Value* pyVec = (lanes > 1) ? emSimd.broadcastF32(pyScalar) : nullptr;
            Value* rowBase = b.build_mul(b.build_add(rowZ, iy), dimX);

            if (lanes > 1) {
                kb.for_range(zero, vecLimit, lanesI32, [&](Value* ix) {
                    Value* ixVec = b.build_vadd(emSimd.broadcastF32(kb.to_f32(ix)), laneRamp);
                    Value* pxVec = emSimd.fma(ixVec, vStepX, vMinX);
                    Value* sdfVec = emitSdfGraphSingleSource(emSimd, graph, root, pxVec, pyVec, pzVec);
                    emSimd.storeLaneF32(outField, b.build_add(rowBase, ix), sdfVec);
                });
            }

            Value* remStart = (lanes > 1) ? b.build_mul(b.build_sdiv(dimX, lanesI32), lanesI32) : zero;
            kb.for_range(remStart, dimX, one, [&](Value* ix) {
                Value* px = b.build_fma(kb.to_f32(ix), stepX, minX);
                Value* sdfVal = emitSdfGraphSingleSource(emScalar, graph, root, px, pyScalar, pzScalar);
                emScalar.storeLaneF32(outField, b.build_add(rowBase, ix), sdfVal);
            });
        });
    });

    b.build_ret_void();
    fn->rebuild_cfg_predecessors();
    return fn;
}

} // namespace mir_detail
#endif // BROMESH_HAS_BRASS_JIT

MirSdfGridKernel compileMirSdfGridKernel(const SdfGraph& graph, int rootNode) {
    MirSdfGridKernel result;
#if BROMESH_HAS_BRASS_JIT
    const int root = (rootNode >= 0) ? rootNode : graph.root();
    const uint32_t lanes = mir_detail::detectHostSimdLanes();
    result.simdLanes = lanes;

    if (Target::host().is_x64() || Target::host().is_aarch64()) {
        try {
            Module cpuMod("mir_sdf_grid_cpu");
            Function* cpuFn = buildGridSdfCpuSimdFunction(cpuMod, graph, root, lanes);
            KernelJit jit(makeSimdJitOptions(lanes));
            result.cpuFn = jit.compile(*cpuFn);
            result.cpuValid = result.cpuFn.is_valid();
        } catch (const std::exception& e) {
            std::cerr << "compileMirSdfGridKernel CPU JIT failed: " << e.what() << std::endl;
        }
    }

    try {
        Module gpuMod("mir_sdf_grid_gpu");
        Function* gpuFn = buildGridSdfGpuFunction(gpuMod, graph, root);
        lowerPortableGpuTargets(gpuMod, *gpuFn, result.gpu, result.spirvKernel);
    } catch (const std::exception& e) {
        result.gpu.portabilityOk = false;
        result.gpu.portabilityDiagnostics += std::string("GPU lowering exception: ") + e.what() + "\n";
    }
#else
    (void)graph;
    (void)rootNode;
#endif
    return result;
}

MirVertexTransformKernel compileMirVertexTransformKernel() {
    MirVertexTransformKernel result;
#if BROMESH_HAS_BRASS_JIT
    const uint32_t lanes = mir_detail::detectHostSimdLanes();
    result.simdLanes = lanes;

    if (Target::host().is_x64() || Target::host().is_aarch64()) {
        try {
            Module cpuMod("mir_vert_xform_cpu");
            Function* cpuFn = buildVertexTransformCpuSimdFunction(cpuMod, lanes);
            KernelJit jit(makeSimdJitOptions(lanes));
            result.cpuFn = jit.compile(*cpuFn);
            result.cpuValid = result.cpuFn.is_valid();
        } catch (const std::exception& e) {
            std::cerr << "compileMirVertexTransformKernel CPU JIT failed: " << e.what() << std::endl;
        }
    }

    try {
        Module gpuMod("mir_vert_xform_gpu");
        Function* gpuFn = buildVertexTransformGpuFunction(gpuMod);
        lowerPortableGpuTargets(gpuMod, *gpuFn, result.gpu, result.spirvKernel);
    } catch (const std::exception& e) {
        result.gpu.portabilityOk = false;
        result.gpu.portabilityDiagnostics += std::string("GPU lowering exception: ") + e.what() + "\n";
    }
#endif
    return result;
}

bool MirSdfGridKernel::evaluateCpuSimd(int dimX, int dimY, int dimZ, const bromath::AABB3& bounds,
                                       std::vector<float>& outField) const {
    if (dimX <= 0 || dimY <= 0 || dimZ <= 0) {
        outField.clear();
        return false;
    }
    const size_t total = static_cast<size_t>(dimX) * static_cast<size_t>(dimY) * static_cast<size_t>(dimZ);
    outField.resize(total, 0.0f);
#if BROMESH_HAS_BRASS_JIT
    if (!cpuValid || !cpuFn.is_valid()) {
        return false;
    }
    float origin[3]{};
    float step[3]{};
    computeGridOriginAndStep(dimX, dimY, dimZ, bounds, origin, step);
    auto fn = cpuFn.as<GridSimdCpuFn>();
    fn(outField.data(), dimX, dimY, dimZ, origin[0], origin[1], origin[2], step[0], step[1], step[2]);
    return true;
#else
    (void)bounds;
    return false;
#endif
}

SDFVolume MirSdfGridKernel::evaluateVolumeCpuSimd(int dimX, int dimY, int dimZ,
                                                  const bromath::AABB3& bounds) const {
    SDFVolume vol;
    if (dimX <= 0 || dimY <= 0 || dimZ <= 0) {
        return vol;
    }
    vol.dimX = dimX;
    vol.dimY = dimY;
    vol.dimZ = dimZ;
    vol.bounds = bounds;
    float origin[3]{};
    computeGridOriginAndStep(dimX, dimY, dimZ, bounds, origin, vol.cellSize);
    evaluateCpuSimd(dimX, dimY, dimZ, bounds, vol.field);
    return vol;
}

bool MirSdfGridKernel::executeVulkan(int dimX, int dimY, int dimZ, const bromath::AABB3& bounds,
                                     std::vector<float>& outField, std::string* errorOut) const {
    if (dimX <= 0 || dimY <= 0 || dimZ <= 0) {
        if (errorOut) *errorOut = "Invalid grid dimensions";
        return false;
    }
#if BROMESH_HAS_BRASS_JIT
    if (!gpu.portabilityOk || spirvKernel.words.empty()) {
        if (errorOut) *errorOut = "SPIR-V kernel is unavailable";
        return false;
    }
    if (!brass::gpu::vulkan_available()) {
        if (errorOut) *errorOut = brass::gpu::vulkan_last_error();
        return false;
    }

    const size_t total = static_cast<size_t>(dimX) * static_cast<size_t>(dimY) * static_cast<size_t>(dimZ);
    const size_t bytes = total * sizeof(float);

    brass::gpu::VulkanBuffer outBuf = brass::gpu::VulkanBuffer::alloc(bytes);
    if (!outBuf.valid()) {
        if (errorOut) *errorOut = brass::gpu::vulkan_last_error();
        return false;
    }

    std::string err;
    brass::gpu::VulkanModule mod = brass::gpu::VulkanModule::load(spirvKernel, &err);
    if (!mod.valid()) {
        if (errorOut) *errorOut = err;
        return false;
    }

    float origin[3]{};
    float step[3]{};
    computeGridOriginAndStep(dimX, dimY, dimZ, bounds, origin, step);

    const uint32_t block = kGpuWorkgroupSize;
    const uint32_t grid = static_cast<uint32_t>((total + block - 1) / block);
    if (!mod.launch_1d(grid, block, {
            outBuf,
            static_cast<int32_t>(dimX),
            static_cast<int32_t>(dimY),
            static_cast<int32_t>(dimZ),
            origin[0], origin[1], origin[2],
            step[0], step[1], step[2]
        }, &err)) {
        if (errorOut) *errorOut = err;
        return false;
    }

    outField.resize(total);
    if (!outBuf.download(outField.data(), bytes)) {
        if (errorOut) *errorOut = brass::gpu::vulkan_last_error();
        return false;
    }
    return true;
#else
    (void)bounds;
    (void)outField;
    if (errorOut) *errorOut = "Brass JIT is disabled";
    return false;
#endif
}

bool MirVertexTransformKernel::executeCpuSimdSoa(const float* inSoa, float* outSoa, const float matrices25[25],
                                                 int32_t numVerts, int32_t strideElems) const {
    if (!inSoa || !outSoa || !matrices25 || numVerts < 0 || strideElems < numVerts) {
        return false;
    }
    if (numVerts == 0) {
        return true;
    }
#if BROMESH_HAS_BRASS_JIT
    if (!cpuValid || !cpuFn.is_valid()) {
        return false;
    }
    auto fn = cpuFn.as<VertTransformCpuFn>();
    fn(outSoa, inSoa, matrices25, numVerts, strideElems);
    return true;
#else
    return false;
#endif
}

bool MirVertexTransformKernel::executeVulkanSoa(const float* inSoa, float* outSoa, const float matrices25[25],
                                                int32_t numVerts, int32_t strideElems,
                                                std::string* errorOut) const {
    if (!inSoa || !outSoa || !matrices25 || numVerts < 0 || strideElems < numVerts) {
        if (errorOut) *errorOut = "Invalid vertex transform arguments";
        return false;
    }
    if (numVerts == 0) {
        return true;
    }
#if BROMESH_HAS_BRASS_JIT
    if (!gpu.portabilityOk || spirvKernel.words.empty()) {
        if (errorOut) *errorOut = "SPIR-V kernel is unavailable";
        return false;
    }
    if (!brass::gpu::vulkan_available()) {
        if (errorOut) *errorOut = brass::gpu::vulkan_last_error();
        return false;
    }

    const size_t soaBytes = static_cast<size_t>(strideElems) * 6u * sizeof(float);
    const size_t matBytes = 25u * sizeof(float);

    brass::gpu::VulkanBuffer inBuf = brass::gpu::VulkanBuffer::alloc(soaBytes);
    brass::gpu::VulkanBuffer outBuf = brass::gpu::VulkanBuffer::alloc(soaBytes);
    brass::gpu::VulkanBuffer matBuf = brass::gpu::VulkanBuffer::alloc(matBytes);
    if (!inBuf.valid() || !outBuf.valid() || !matBuf.valid()) {
        if (errorOut) *errorOut = brass::gpu::vulkan_last_error();
        return false;
    }
    if (!inBuf.upload(inSoa, soaBytes) || !matBuf.upload(matrices25, matBytes)) {
        if (errorOut) *errorOut = brass::gpu::vulkan_last_error();
        return false;
    }

    std::string err;
    brass::gpu::VulkanModule mod = brass::gpu::VulkanModule::load(spirvKernel, &err);
    if (!mod.valid()) {
        if (errorOut) *errorOut = err;
        return false;
    }

    const uint32_t block = kGpuWorkgroupSize;
    const uint32_t grid = static_cast<uint32_t>((static_cast<size_t>(numVerts) + block - 1) / block);
    if (!mod.launch_1d(grid, block, {outBuf, inBuf, matBuf, numVerts, strideElems}, &err)) {
        if (errorOut) *errorOut = err;
        return false;
    }

    if (!outBuf.download(outSoa, soaBytes)) {
        if (errorOut) *errorOut = brass::gpu::vulkan_last_error();
        return false;
    }
    return true;
#else
    if (errorOut) *errorOut = "Brass JIT is disabled";
    return false;
#endif
}

bool MirVertexTransformKernel::transformMeshCpuSimd(const MeshData& src, MeshData& dst, const float mat4[16]) const {
    return transformMeshViaSoa(src, dst, mat4, [&](const float* inSoa, float* outSoa, const float* mats,
                                                   int32_t numVerts, int32_t stride) {
        return executeCpuSimdSoa(inSoa, outSoa, mats, numVerts, stride);
    });
}

bool MirVertexTransformKernel::transformMeshVulkan(const MeshData& src, MeshData& dst, const float mat4[16],
                                                   std::string* errorOut) const {
    return transformMeshViaSoa(src, dst, mat4, [&](const float* inSoa, float* outSoa, const float* mats,
                                                   int32_t numVerts, int32_t stride) {
        return executeVulkanSoa(inSoa, outSoa, mats, numVerts, stride, errorOut);
    });
}

bool isMirVulkanAvailable() {
#if BROMESH_HAS_BRASS_JIT
    return brass::gpu::vulkan_available();
#else
    return false;
#endif
}

std::string mirVulkanDeviceSummary() {
#if BROMESH_HAS_BRASS_JIT
    if (!brass::gpu::vulkan_available()) {
        return "unavailable: " + brass::gpu::vulkan_last_error();
    }
    const brass::gpu::VulkanDeviceCaps* caps = brass::gpu::vulkan_device_caps();
    if (!caps) {
        return "available (no device caps)";
    }
    return caps->info.name + " (" + caps->info.type + ", driver " + caps->driver_name + ")";
#else
    return "unavailable (Brass JIT disabled)";
#endif
}

} // namespace bromesh
