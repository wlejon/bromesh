#include "bromesh/isosurface/jit/mir_mesh_kernels.h"
#include "bromesh/isosurface/jit/sdf_compiler.h"
#include "bromesh/isosurface/jit/jit_mesher.h"
#include "bromesh/manipulation/normals.h"
#include "bromesh/manipulation/transform.h"
#include "bromesh/primitives/primitives.h"

#if BROMESH_HAS_BRASS_JIT
#include <brass/target/spirv/spirv_binary_reader.hpp>
#include <brass/verify/gpu_kernel_prover.hpp>
#endif

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

namespace {

int gChecksRun = 0;
int gChecksPassed = 0;

#define CHECK(cond, msg) do { \
    ++gChecksRun; \
    if (!(cond)) { \
        std::fprintf(stderr, "FAIL: %s (line %d)\n", msg, __LINE__); \
    } else { \
        ++gChecksPassed; \
    } \
} while (0)

void verifyGpuArtifacts(const bromesh::MirGpuArtifacts& gpu, const char* expectedEntry) {
    CHECK(gpu.portabilityOk, "GPU MIR must pass brass::verify_gpu_portability");
    if (!gpu.portabilityOk) {
        std::fprintf(stderr, "  portability diagnostics: %s\n", gpu.portabilityDiagnostics.c_str());
    }

    CHECK(!gpu.spirvWords.empty(), "SPIR-V output must be non-empty");
    CHECK(gpu.spirvWords.size() >= 5 && gpu.spirvWords[0] == 0x07230203u,
          "SPIR-V binary must start with magic 0x07230203");
    CHECK(gpu.spirvEntry == expectedEntry, "SPIR-V entry name must match");

#if BROMESH_HAS_BRASS_JIT
    std::string spvErr;
    CHECK(brass::spirv::verify_spirv_module(gpu.spirvWords, &spvErr),
          "SPIR-V module must pass verify_spirv_module");
    if (!spvErr.empty()) {
        std::fprintf(stderr, "  spirv verify err: %s\n", spvErr.c_str());
    }

    brass::spirv::Module spvMod;
    std::string parseErr;
    CHECK(brass::spirv::parse_binary(gpu.spirvWords, spvMod, &parseErr),
          "SPIR-V binary must parse for GPU safety audit");
    auto audit = brass::verify::prove_spirv_shader_safety(spvMod);
    CHECK(audit.structural_ok && audit.dominance_ok && audit.barriers_uniform &&
          audit.shared_bounds_ok && audit.div_sqrt_ok && audit.hazards.empty(),
          "SPIR-V shader must pass GPU safety prover");
#endif

    CHECK(!gpu.ptxSource.empty(), "PTX source must be non-empty");
    const std::string ptxEntry = std::string(".visible .entry ") + expectedEntry;
    CHECK(gpu.ptxSource.find(ptxEntry) != std::string::npos,
          "PTX must expose expected .visible .entry");
    CHECK(gpu.ptxSource.find("st.global") != std::string::npos,
          "PTX must emit global memory store");

    CHECK(!gpu.mslSource.empty(), "MSL source must be non-empty");
    CHECK(gpu.mslSource.find("#include <metal_stdlib>") != std::string::npos,
          "MSL must include <metal_stdlib>");
    CHECK(gpu.mslSource.find(expectedEntry) != std::string::npos,
          "MSL must contain expected kernel entry name");
}

float maxAbsDiff(const std::vector<float>& a, const std::vector<float>& b) {
    if (a.size() != b.size()) return 1e9f;
    float maxD = 0.0f;
    for (size_t i = 0; i < a.size(); ++i) {
        const float d = std::fabs(a[i] - b[i]);
        if (d > maxD || !std::isfinite(d)) {
            maxD = std::isfinite(d) ? d : 1e9f;
        }
    }
    return maxD;
}

void testGridSdfCsgKernel() {
    bromesh::SdfGraph g;
    int sp = g.translate(g.sphere(0.85f), {-0.35f, 0.10f, 0.0f});
    int bx = g.rotateY(g.roundedBox({0.55f, 0.45f, 0.50f}, 0.08f), 0.40f);
    int u1 = g.opSmoothUnion(sp, bx, 0.20f);
    int tr = g.translate(g.torus(0.75f, 0.18f), {0.25f, -0.15f, 0.10f});
    int u2 = g.opSmoothUnion(u1, tr, 0.15f);
    int cy = g.cylinder(0.25f, 1.10f);
    g.opSmoothSubtraction(u2, cy, 0.08f);

    bromesh::MirSdfGridKernel kernel = bromesh::compileMirSdfGridKernel(g);
    CHECK(kernel.cpuValid, "Grid SDF CPU SIMD kernel must compile");
    CHECK(kernel.simdLanes >= 4, "Host CPU SIMD kernel must use >= 4 SIMD lanes");
    verifyGpuArtifacts(kernel.gpu, "eval_sdf_grid_gpu");

    // Non-multiple-of-8 dimensions exercise both SIMD vector body and scalar tail.
    const int dimX = 29;
    const int dimY = 27;
    const int dimZ = 25;
    const bromath::AABB3 bounds{{-1.6f, -1.6f, -1.6f}, {1.6f, 1.6f, 1.6f}};

    std::vector<float> refField;
    std::vector<float> cpuField;
    bromesh::evaluateSdfGridReference(g, dimX, dimY, dimZ, bounds, refField);
    CHECK(kernel.evaluateCpuSimd(dimX, dimY, dimZ, bounds, cpuField),
          "Grid SDF CPU SIMD evaluation must succeed");
    CHECK(cpuField.size() == refField.size(), "Grid SDF CPU SIMD field size must match reference");

    const float cpuMaxDiff = maxAbsDiff(cpuField, refField);
    std::printf("[GRID SDF CSG] CPU SIMD (lanes=%u) vs C++ ref maxDiff = %.6g\n",
                kernel.simdLanes, cpuMaxDiff);
    CHECK(cpuMaxDiff < 1e-5f, "CPU SIMD grid SDF must match C++ reference within 1e-5");

    if (bromesh::isMirVulkanAvailable()) {
        std::vector<float> vkField;
        std::string vkErr;
        const bool ok = kernel.executeVulkan(dimX, dimY, dimZ, bounds, vkField, &vkErr);
        if (!ok) {
            std::fprintf(stderr, "  Vulkan grid SDF error: %s\n", vkErr.c_str());
        }
        CHECK(ok, "Vulkan grid SDF execution must succeed when Vulkan is available");
        const float vkMaxDiff = maxAbsDiff(vkField, cpuField);
        std::printf("[GRID SDF CSG] Vulkan (%s) vs CPU SIMD maxDiff = %.6g\n",
                    bromesh::mirVulkanDeviceSummary().c_str(), vkMaxDiff);
        CHECK(vkMaxDiff < 1e-4f, "Vulkan grid SDF must match CPU SIMD within 1e-4");
    } else {
        std::printf("[SKIP] Vulkan unavailable (%s)\n", bromesh::mirVulkanDeviceSummary().c_str());
    }
}

void testGridSdfNoiseAndTrigKernel() {
    bromesh::SdfGraph g;
    int sp = g.sphere(0.95f);
    int noise = g.noise3D({2.2f, 2.2f, 2.2f}, 0.06f);
    int disp = g.displace(sp, noise);
    int wave = g.customExpr("sin(x * 3.0) * cos(y * 3.0) * 0.04", 0.5f);
    g.displace(disp, wave);

    bromesh::MirSdfGridKernel kernel = bromesh::compileMirSdfGridKernel(g);
    CHECK(kernel.cpuValid, "Noise+Trig SDF CPU SIMD kernel must compile");
    verifyGpuArtifacts(kernel.gpu, "eval_sdf_grid_gpu");

    const int dimX = 21;
    const int dimY = 19;
    const int dimZ = 17;
    const bromath::AABB3 bounds{{-1.4f, -1.4f, -1.4f}, {1.4f, 1.4f, 1.4f}};

    std::vector<float> refField;
    std::vector<float> cpuField;
    bromesh::evaluateSdfGridReference(g, dimX, dimY, dimZ, bounds, refField);
    CHECK(kernel.evaluateCpuSimd(dimX, dimY, dimZ, bounds, cpuField),
          "Noise+Trig CPU SIMD evaluation must succeed");

    const float cpuMaxDiff = maxAbsDiff(cpuField, refField);
    std::printf("[GRID SDF NOISE+TRIG] CPU SIMD (lanes=%u) vs C++ ref maxDiff = %.6g\n",
                kernel.simdLanes, cpuMaxDiff);
    CHECK(cpuMaxDiff < 1e-5f, "Noise+Trig CPU SIMD SDF must match C++ reference within 1e-5");

    if (bromesh::isMirVulkanAvailable()) {
        std::vector<float> vkField;
        std::string vkErr;
        const bool ok = kernel.executeVulkan(dimX, dimY, dimZ, bounds, vkField, &vkErr);
        if (!ok) {
            std::fprintf(stderr, "  Vulkan noise+trig SDF error: %s\n", vkErr.c_str());
        }
        CHECK(ok, "Vulkan noise+trig grid SDF execution must succeed");
        const float vkMaxDiff = maxAbsDiff(vkField, cpuField);
        std::printf("[GRID SDF NOISE+TRIG] Vulkan vs CPU SIMD maxDiff = %.6g\n", vkMaxDiff);
        CHECK(vkMaxDiff < 1e-4f, "Vulkan noise+trig SDF must match CPU SIMD within 1e-4");
    }
}

void testVertexTransformAndNormalKernel() {
    bromesh::MeshData src = bromesh::sphere(1.0f, 16, 12);
    // Append a couple of custom vertices (non-multiple-of-8 total + zero-normal edge case).
    src.positions.push_back(1.25f); src.positions.push_back(-0.5f); src.positions.push_back(0.75f);
    src.normals.push_back(0.0f); src.normals.push_back(0.0f); src.normals.push_back(0.0f);
    src.positions.push_back(-0.3f); src.positions.push_back(0.9f); src.positions.push_back(-1.1f);
    src.normals.push_back(3.0f); src.normals.push_back(4.0f); src.normals.push_back(0.0f);

    // Non-uniform scale + rotation around Z + translation (column-major 4x4)
    const float theta = 0.55f;
    const float c = std::cos(theta);
    const float s = std::sin(theta);
    const float sx = 1.50f, sy = 0.75f, sz = 1.25f;
    const float mat4[16] = {
        c * sx,  s * sx, 0.0f, 0.0f,
       -s * sy,  c * sy, 0.0f, 0.0f,
        0.0f,    0.0f,   sz,   0.0f,
        0.40f,  -0.25f,  0.80f, 1.0f
    };

    bromesh::MirVertexTransformKernel kernel = bromesh::compileMirVertexTransformKernel();
    CHECK(kernel.cpuValid, "Vertex transform CPU SIMD kernel must compile");
    verifyGpuArtifacts(kernel.gpu, "transform_mesh_gpu");

    bromesh::MeshData refMesh;
    bromesh::MeshData cpuMesh;
    bromesh::transformMeshReference(src, refMesh, mat4);
    CHECK(kernel.transformMeshCpuSimd(src, cpuMesh, mat4),
          "Vertex transform CPU SIMD execution must succeed");

    const float posDiff = maxAbsDiff(cpuMesh.positions, refMesh.positions);
    const float nrmDiff = maxAbsDiff(cpuMesh.normals, refMesh.normals);
    std::printf("[VERTEX XFORM] verts=%zu CPU SIMD (lanes=%u) vs C++ ref: posDiff=%.6g normDiff=%.6g\n",
                src.positions.size() / 3, kernel.simdLanes, posDiff, nrmDiff);
    CHECK(posDiff < 1e-5f, "CPU SIMD vertex positions must match C++ reference within 1e-5");
    CHECK(nrmDiff < 1e-5f, "CPU SIMD normalized normals must match C++ reference within 1e-5");

    if (bromesh::isMirVulkanAvailable()) {
        bromesh::MeshData vkMesh;
        std::string vkErr;
        const bool ok = kernel.transformMeshVulkan(src, vkMesh, mat4, &vkErr);
        if (!ok) {
            std::fprintf(stderr, "  Vulkan vertex transform error: %s\n", vkErr.c_str());
        }
        CHECK(ok, "Vulkan vertex transform execution must succeed");
        const float vkPosDiff = maxAbsDiff(vkMesh.positions, cpuMesh.positions);
        const float vkNrmDiff = maxAbsDiff(vkMesh.normals, cpuMesh.normals);
        std::printf("[VERTEX XFORM] Vulkan vs CPU SIMD: posDiff=%.6g normDiff=%.6g\n",
                    vkPosDiff, vkNrmDiff);
        CHECK(vkPosDiff < 1e-4f, "Vulkan vertex positions must match CPU SIMD within 1e-4");
        CHECK(vkNrmDiff < 1e-4f, "Vulkan normalized normals must match CPU SIMD within 1e-4");
    }
}

void testTransformMeshSimdWiring() {
    bromesh::MeshData base = bromesh::sphere(1.0f, 16, 12);
    bromesh::generateTangents(base);
    CHECK(base.vertexCount() >= 64 && base.hasNormals() && base.hasTangents(),
          "base sphere for transformMesh test has >= 64 vertices, normals, and tangents");

    const float posMat[16] = {
        1.2f,  0.3f, 0.0f, 0.0f,
       -0.2f,  0.9f, 0.1f, 0.0f,
        0.0f, -0.1f, 1.4f, 0.0f,
        0.5f, -0.3f, 1.1f, 1.0f
    };
    const float negMat[16] = {
       -1.1f,  0.2f, 0.0f, 0.0f,
        0.1f,  0.8f, 0.0f, 0.0f,
        0.0f,  0.0f, 1.3f, 0.0f,
       -0.4f,  0.7f, 0.2f, 1.0f
    };

    auto verifyTransformCase = [&](const bromesh::MeshData& src, const float m[16], const char* label) {
        bromesh::MeshData actual = src;
        bromesh::transformMesh(actual, m);

        bromesh::MeshData ref;
        bromesh::transformMeshReference(src, ref, m);

        const float pDiff = maxAbsDiff(actual.positions, ref.positions);
        CHECK(pDiff < 1e-5f, label);
        if (src.hasNormals()) {
            const float nDiff = maxAbsDiff(actual.normals, ref.normals);
            CHECK(nDiff < 1e-5f, label);
        } else {
            CHECK(actual.normals.empty(), "transformMesh must not allocate normals when absent");
        }
        if (src.hasTangents()) {
            const float det =
                m[0] * (m[5] * m[10] - m[6] * m[9]) -
                m[4] * (m[1] * m[10] - m[2] * m[9]) +
                m[8] * (m[1] * m[6] - m[2] * m[5]);
            const float wantW = (det < 0.0f ? -1.0f : 1.0f) * src.tangents[3];
            CHECK(std::fabs(actual.tangents[3] - wantW) < 1e-5f,
                  "transformMesh tangent handedness must track det < 0");
        }
    };

    verifyTransformCase(base, posMat, "transformMesh positive-det matches ref");
    verifyTransformCase(base, negMat, "transformMesh negative-det mirror matches ref");

    bromesh::MeshData noNormals = base;
    noNormals.normals.clear();
    noNormals.tangents.clear();
    verifyTransformCase(noNormals, posMat, "transformMesh without normals matches ref");
}

void testSdfCompilerSimdMslBundlesAndMesher() {
    bromesh::SdfGraph g;
    const int sph = g.sphere(0.9f);
    const int bx = g.translate(g.box({0.55f, 0.55f, 0.55f}), {0.65f, 0.0f, 0.0f});
    g.setRoot(g.opSmoothUnion(sph, bx, 0.2f));

    auto& comp = bromesh::JitSdfCompiler::instance();
    const bromesh::CompiledSdfKernel ck = comp.compile(g);
    CHECK(ck.isJit, "JitSdfCompiler::compile must produce valid JIT kernels");
    CHECK(ck.simdGridFn != nullptr, "JitSdfCompiler::compile must expose eval_grid_simd");

    const bromath::AABB3 bounds{{-2.0f, -2.0f, -2.0f}, {2.0f, 2.0f, 2.0f}};
    const int dim = 24;
    const bromesh::SDFVolume unprunedVol = comp.evaluateVolume(g, dim, dim, dim, bounds, -1, true, 0.0f, false);
    const bromesh::SDFVolume prunedVol = comp.evaluateVolume(g, dim, dim, dim, bounds, -1, true, 0.0f, true);
    std::vector<float> refField;
    bromesh::evaluateSdfGridReference(g, dim, dim, dim, bounds, refField);

    const float unprunedDiff = maxAbsDiff(unprunedVol.field, refField);
    CHECK(unprunedDiff < 1e-5f, "evaluateVolume (unpruned SIMD) must match scalar ref");

    const bromesh::MeshData mcField = bromesh::marchingCubesFromField(
        g, unprunedVol.field.data(), dim, dim, dim, bounds, 0.0f, true, true, nullptr);
    const bromesh::MeshData mcSdf = bromesh::marchingCubesFromSDF(
        g, dim, dim, dim, bounds, 0.0f, true, true, false);
    CHECK(!mcField.empty() && mcField.positions.size() == mcSdf.positions.size(),
          "marchingCubesFromField matches marchingCubesFromSDF vertex count");
    CHECK(maxAbsDiff(mcField.positions, mcSdf.positions) < 1e-5f,
          "marchingCubesFromField matches marchingCubesFromSDF positions");
    CHECK(!prunedVol.field.empty(), "evaluateVolume (pruned) succeeds");

    const bromesh::MeshData snField = bromesh::surfaceNetsFromField(
        g, unprunedVol.field.data(), dim, dim, dim, bounds, 0.0f, nullptr);
    const bromesh::MeshData snSdf = bromesh::surfaceNetsFromSDF(
        g, dim, dim, dim, bounds, 0.0f, false);
    CHECK(!snField.empty() && snField.positions.size() == snSdf.positions.size(),
          "surfaceNetsFromField matches surfaceNetsFromSDF vertex count");
    CHECK(maxAbsDiff(snField.positions, snSdf.positions) < 1e-5f,
          "surfaceNetsFromField matches surfaceNetsFromSDF positions");

    // Multi-target SPIR-V + MSL + bundle round-trip
    const auto spv1 = comp.compile_spirv(g);
    const auto spv2 = comp.compile_spirv(g);
    CHECK(!spv1.empty() && spv1 == spv2, "compile_spirv hit matches initial compilation");

    const std::string msl = comp.compile_msl(g);
    CHECK(!msl.empty() && msl.find("#include <metal_stdlib>") != std::string::npos &&
              msl.find("eval_sdf_gpu") != std::string::npos,
          "compile_msl emits valid Metal compute shader");

    const std::filesystem::path bundlePath =
        std::filesystem::temp_directory_path() / "bromesh_test_sdf_bundle.bacb";
    std::filesystem::remove(bundlePath);
    std::string bundleErr;
    CHECK(comp.export_kernel_bundle(bundlePath.string(), &bundleErr),
          "export_kernel_bundle must succeed");
    CHECK(std::filesystem::exists(bundlePath) && std::filesystem::file_size(bundlePath) > 0,
          "exported bundle file must be non-empty");
    comp.clearCache();
    CHECK(comp.import_kernel_bundle(bundlePath.string(), &bundleErr),
          "import_kernel_bundle must succeed");
    const auto spvImported = comp.compile_spirv(g);
    CHECK(spvImported == spv1, "re-imported SPIR-V from bundle matches original");
    std::filesystem::remove(bundlePath);
}

} // namespace

int main() {
    testGridSdfCsgKernel();
    testGridSdfNoiseAndTrigKernel();
    testVertexTransformAndNormalKernel();
    testTransformMeshSimdWiring();
    testSdfCompilerSimdMslBundlesAndMesher();
    std::printf("mir_mesh_kernels tests: %d/%d passed\n", gChecksPassed, gChecksRun);
    return (gChecksPassed == gChecksRun) ? 0 : 1;
}
