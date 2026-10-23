#pragma once

#include "bromesh/isosurface/jit/sdf_node.h"
#include "bromesh/analysis/sdf.h"
#include <bromath/aabb.h>
#include <bromath/vec.h>

#if BROMESH_HAS_BRASS_JIT
#include <brass/codegen/kernel_jit.hpp>
#include <brass/codegen/sdf_builder.hpp>
#include <brass/target/spirv_target.hpp>
#include <brass/target/ptx_target.hpp>
#endif

#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace bromesh {

typedef void (*SdfGridEvalFn)(float* out_field, const float* bounds_min, const float* cell_size, int32_t dim_x, int32_t dim_y, int32_t dim_z);
typedef void (*SdfSpanEvalFn)(float* out_row, float ox, float py, float pz, float sx, int64_t x_start, int64_t x_end);
typedef float (*SdfPointEvalFn)(float x, float y, float z);
typedef float (*SdfGradEvalFn)(float x, float y, float z, float* out_gx, float* out_gy, float* out_gz);
typedef void (*SdfBlockBoundsFn)(float min_x, float min_y, float min_z, float max_x, float max_y, float max_z, float* out_dmin, float* out_dmax);

struct SdfContact {
    float distance = 0.0f;
    bromath::Vec3 normal{0.0f, 1.0f, 0.0f};
    bromath::Vec3 closestSurfacePoint{0.0f, 0.0f, 0.0f};
    bromath::Vec3 gradient{0.0f, 0.0f, 0.0f};
};

struct SdfInterval {
    float lo = 0.0f;
    float hi = 0.0f;
};

/// Compute guaranteed [dMin, dMax] interval / Lipschitz bounds over blockBounds for nodeId.
/// Supports CSG branch pruning: when one branch of a CSG operator provably dominates the other
/// across blockBounds, *activeNodeOut (if non-null) receives the dominating subtree node ID.
void evalBlockBounds(const SdfGraph& graph, int nodeId, const bromath::AABB3& blockBounds,
                     float& dMin, float& dMax, int* activeNodeOut = nullptr);

/// Return a pruned copy of graph for a macro-block where dominated CSG branches are eliminated.
SdfGraph pruneCsgGraphForBlock(const SdfGraph& graph, const bromath::AABB3& blockBounds, int rootNode = -1);

struct CompiledSdfKernel {
#if BROMESH_HAS_BRASS_JIT
    brass::codegen::KernelFunction gridKernel;
    brass::codegen::KernelFunction spanKernel;
    brass::codegen::KernelFunction pointKernel;
    brass::codegen::KernelFunction gradKernel;
    brass::codegen::KernelFunction blockBoundsKernel;
#endif
    SdfGridEvalFn gridFn = nullptr;
    SdfSpanEvalFn spanFn = nullptr;
    SdfPointEvalFn pointFn = nullptr;
    SdfGradEvalFn gradFn = nullptr;
    SdfBlockBoundsFn blockBoundsFn = nullptr;
    uint64_t hash = 0;
    bool isJit = false;

    /// O(1) SDF contact query at (x, y, z) using 1-pass forward-mode symbolic auto-diff.
    SdfContact queryContact(float x, float y, float z) const;

    /// Evaluate [dMin, dMax] interval bounds over a 3D box using compiled interval MIR.
    void evalBlockBounds(const bromath::AABB3& blockBounds, float& dMin, float& dMax) const;
};

class JitSdfCompiler {
public:
    JitSdfCompiler() = default;
    ~JitSdfCompiler() = default;

    static JitSdfCompiler& instance();

    CompiledSdfKernel compile(const SdfGraph& graph, int rootNode = -1);

    std::vector<uint32_t> compile_spirv(const SdfGraph& graph, int rootNode = -1);
    std::string compile_ptx(const SdfGraph& graph, int rootNode = -1);

    std::vector<uint32_t> compileSpirv(const SdfGraph& graph, int rootNode = -1) {
        return compile_spirv(graph, rootNode);
    }
    std::string compilePtx(const SdfGraph& graph, int rootNode = -1) {
        return compile_ptx(graph, rootNode);
    }

    SDFVolume evaluateVolume(const SdfGraph& graph, int dimX, int dimY, int dimZ,
                             const bromath::AABB3& bounds, int rootNode = -1, bool parallel = true,
                             float isoLevel = 0.0f, bool enableMacroBlockPruning = true);

    void evalBlockBounds(const SdfGraph& graph, int nodeId, const bromath::AABB3& blockBounds,
                         float& dMin, float& dMax, int* activeNodeOut = nullptr) const {
        ::bromesh::evalBlockBounds(graph, nodeId, blockBounds, dMin, dMax, activeNodeOut);
    }

    void clearCache();

private:
    std::mutex m_cacheMutex;
    std::unordered_map<uint64_t, CompiledSdfKernel> m_cache;
};

using SdfCompiler = JitSdfCompiler;

} // namespace bromesh
