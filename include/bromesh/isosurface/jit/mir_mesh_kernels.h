#pragma once

#include "bromesh/analysis/sdf.h"
#include "bromesh/isosurface/jit/sdf_node.h"
#include "bromesh/mesh_data.h"
#include <bromath/aabb.h>

#if BROMESH_HAS_BRASS_JIT
#include <brass/codegen/kernel_jit.hpp>
#include <brass/target/msl_target.hpp>
#include <brass/target/ptx_target.hpp>
#include <brass/target/spirv_target.hpp>
#endif

#include <cstdint>
#include <string>
#include <vector>

namespace bromesh {

/// Artifacts produced by compiling a single-source GPU Brass MIR kernel across
/// SPIR-V (Vulkan), PTX (CUDA), and MSL (Metal), after verify_gpu_portability.
struct MirGpuArtifacts {
    bool portabilityOk = false;
    std::string portabilityDiagnostics;
    std::vector<uint32_t> spirvWords;
    std::string spirvEntry;
    std::string ptxSource;
    std::string mslSource;
    std::string mslEntry;
    uint32_t workgroupX = 64;
};

/// Single-source Brass MIR kernel for batched 3D regular-grid SDF evaluation.
/// Lowers an SdfGraph into:
///   - CPU SIMD MIR (AVX2 f32x8 on x86-64, f32x4 on SSE/AArch64) with scalar tail
///   - Portable GPU MIR verified via brass::verify_gpu_portability and compiled
///     to SPIR-V, PTX, and Metal MSL (executable via VulkanModule when available).
struct MirSdfGridKernel {
    uint32_t simdLanes = 1;
    bool cpuValid = false;
    MirGpuArtifacts gpu;

#if BROMESH_HAS_BRASS_JIT
    brass::codegen::KernelFunction cpuFn;
    brass::target::SpirvKernel spirvKernel;
#endif

    bool evaluateCpuSimd(int dimX, int dimY, int dimZ, const bromath::AABB3& bounds,
                         std::vector<float>& outField) const;

    SDFVolume evaluateVolumeCpuSimd(int dimX, int dimY, int dimZ, const bromath::AABB3& bounds) const;

    bool executeVulkan(int dimX, int dimY, int dimZ, const bromath::AABB3& bounds,
                       std::vector<float>& outField, std::string* errorOut = nullptr) const;
};

/// Single-source Brass MIR kernel for batched SoA vertex affine transform
/// and normal matrix transform + normalization.
/// Buffer layout for inSoa / outSoa: 6 contiguous channels of strideElems floats
///   [0]: px, [1]: py, [2]: pz, [3]: nx, [4]: ny, [5]: nz.
/// Matrix buffer layout (25 floats):
///   [0..15]: column-major 4x4 vertex transform matrix
///   [16..24]: column-major 3x3 normal matrix (inverse-transpose of upper 3x3).
struct MirVertexTransformKernel {
    uint32_t simdLanes = 1;
    bool cpuValid = false;
    MirGpuArtifacts gpu;

#if BROMESH_HAS_BRASS_JIT
    brass::codegen::KernelFunction cpuFn;
    brass::target::SpirvKernel spirvKernel;
#endif

    bool executeCpuSimdSoa(const float* inSoa, float* outSoa, const float matrices25[25],
                           int32_t numVerts, int32_t strideElems) const;

    bool executeVulkanSoa(const float* inSoa, float* outSoa, const float matrices25[25],
                          int32_t numVerts, int32_t strideElems, std::string* errorOut = nullptr) const;

    bool transformMeshCpuSimd(const MeshData& src, MeshData& dst, const float mat4[16]) const;

    bool transformMeshVulkan(const MeshData& src, MeshData& dst, const float mat4[16],
                             std::string* errorOut = nullptr) const;
};

/// Compile a single-source CPU SIMD + portable GPU SDF grid evaluation kernel.
MirSdfGridKernel compileMirSdfGridKernel(const SdfGraph& graph, int rootNode = -1);

/// Compile a single-source CPU SIMD + portable GPU vertex transform & normal normalization kernel.
MirVertexTransformKernel compileMirVertexTransformKernel();

/// Query whether Brass's Vulkan compute runtime is available on this host.
bool isMirVulkanAvailable();

/// Human-readable Vulkan device summary (or unavailable reason).
std::string mirVulkanDeviceSummary();

/// C++ scalar reference implementation for regular-grid SDF evaluation.
void evaluateSdfGridReference(const SdfGraph& graph, int dimX, int dimY, int dimZ,
                              const bromath::AABB3& bounds, std::vector<float>& outField,
                              int rootNode = -1);

/// Compute column-major 3x3 normal matrix (inverse-transpose of upper 3x3 of mat4).
void computeNormalMatrix3x3(const float mat4[16], float outNormMat3[9]);

/// C++ reference vertex affine transform + normal transform & normalization on MeshData.
void transformMeshReference(const MeshData& src, MeshData& dst, const float mat4[16]);

} // namespace bromesh
