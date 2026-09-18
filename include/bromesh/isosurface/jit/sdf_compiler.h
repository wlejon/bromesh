#pragma once

#include "bromesh/isosurface/jit/sdf_node.h"
#include "bromesh/analysis/sdf.h"
#include <bromath/aabb.h>

#if BROMESH_HAS_BRASS_JIT
#include <brass/codegen/kernel_jit.hpp>
#include <brass/codegen/sdf_builder.hpp>
#endif

#include <cstdint>
#include <memory>
#include <mutex>
#include <unordered_map>

namespace bromesh {

typedef void (*SdfGridEvalFn)(float* out_field, const float* bounds_min, const float* cell_size, int32_t dim_x, int32_t dim_y, int32_t dim_z);
typedef float (*SdfPointEvalFn)(float x, float y, float z);

struct CompiledSdfKernel {
#if BROMESH_HAS_BRASS_JIT
    brass::codegen::KernelFunction gridKernel;
    brass::codegen::KernelFunction pointKernel;
#endif
    SdfGridEvalFn gridFn = nullptr;
    SdfPointEvalFn pointFn = nullptr;
    uint64_t hash = 0;
    bool isJit = false;
};

class JitSdfCompiler {
public:
    JitSdfCompiler() = default;
    ~JitSdfCompiler() = default;

    static JitSdfCompiler& instance();

    CompiledSdfKernel compile(const SdfGraph& graph, int rootNode = -1);

    SDFVolume evaluateVolume(const SdfGraph& graph, int dimX, int dimY, int dimZ,
                             const bromath::AABB3& bounds, int rootNode = -1);

    void clearCache();

private:
    std::mutex m_cacheMutex;
    std::unordered_map<uint64_t, CompiledSdfKernel> m_cache;
};

} // namespace bromesh
