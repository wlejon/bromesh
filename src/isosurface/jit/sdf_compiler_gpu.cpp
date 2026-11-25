#include "sdf_compiler_internal.h"

#if BROMESH_HAS_BRASS_JIT
#include <brass/cache/artifact_cache.hpp>
#include <brass/core/diagnostics.hpp>
#include <brass/mir/gpu_portability.hpp>
#include <brass/target/msl_target.hpp>
#include <brass/target/ptx_target.hpp>
#include <brass/target/spirv_target.hpp>
#endif

#include <cstring>
#include <filesystem>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace bromesh {

#if BROMESH_HAS_BRASS_JIT
namespace {

uint64_t gpu_mem_key(const SdfGraph& graph, int rootNode) {
    uint64_t h = graph.computeHash();
    const int root = (rootNode >= 0) ? rootNode : graph.root();
    h ^= (static_cast<uint64_t>(static_cast<uint32_t>(root)) * 0x9e3779b97f4a7c15ULL);
    return h;
}

struct GpuTargetMemCache {
    std::mutex mu;
    std::unordered_map<uint64_t, std::vector<uint32_t>> spirv;
    std::unordered_map<uint64_t, std::string> ptx;
    std::unordered_map<uint64_t, std::string> msl;
};

GpuTargetMemCache& gpu_mem_cache() {
    static auto* c = new GpuTargetMemCache();
    return *c;
}

brass::cache::CacheKey build_gpu_artifact_key(const char* domain, const SdfGraph& graph, int root,
                                            bool loadsThroughVar = false) {
    brass::cache::KeyBuilder kb(domain);
    kb.add_brass_codegen();
    append_sdf_graph_to_key(kb, graph, root);
    kb.add_bool("loads_through_var", loadsThroughVar);
    return kb.finish();
}

} // namespace

void clearGpuSdfCompilerCaches() {
    auto& c = gpu_mem_cache();
    std::lock_guard<std::mutex> lk(c.mu);
    c.spirv.clear();
    c.ptx.clear();
    c.msl.clear();
}
#endif

std::vector<uint32_t> JitSdfCompiler::compile_spirv(const SdfGraph& graph, int rootNode) {
#if BROMESH_HAS_BRASS_JIT
    if (graph.empty()) return {};
    const int root = (rootNode >= 0) ? rootNode : graph.root();
    const uint64_t mkey = gpu_mem_key(graph, root);
    auto& mc = gpu_mem_cache();
    {
        std::lock_guard<std::mutex> lk(mc.mu);
        if (auto it = mc.spirv.find(mkey); it != mc.spirv.end()) {
            return it->second;
        }
    }

#ifdef __APPLE__
    constexpr bool kLoadsThroughVar = true;
#else
    constexpr bool kLoadsThroughVar = false;
#endif
    const brass::cache::CacheKey dkey = build_gpu_artifact_key("bromesh.sdf.spirv/1", graph, root, kLoadsThroughVar);
    auto& store = bromesh_artifact_cache();
    if (auto art = store.get(dkey)) {
        const auto words = art->as<uint32_t>("spirv");
        if (!words.empty()) {
            std::vector<uint32_t> out(words.begin(), words.end());
            std::lock_guard<std::mutex> lk(mc.mu);
            mc.spirv[mkey] = out;
            return out;
        }
    }

    brass::Module mod("sdf_spirv_" + std::to_string(mkey));
    brass::Function* fn = buildGpuSdfKernel(mod, graph, root);
    brass::DiagnosticReporter portDiag;
    (void)brass::verify_gpu_portability(*fn, {}, &portDiag);

    brass::target::SpirvOptions spvOpts;
    spvOpts.loads_through_variable = kLoadsThroughVar;
    auto kernel = brass::target::SpirvTarget::compile(*fn, spvOpts);
    if (!kernel.words.empty()) {
        std::string putErr;
        store.put(dkey, {brass::cache::BlobRef::of("spirv", std::span<const uint32_t>(kernel.words))}, &putErr);
        std::lock_guard<std::mutex> lk(mc.mu);
        mc.spirv[mkey] = kernel.words;
    }
    return kernel.words;
#else
    (void)graph;
    (void)rootNode;
    return {};
#endif
}

std::string JitSdfCompiler::compile_ptx(const SdfGraph& graph, int rootNode) {
#if BROMESH_HAS_BRASS_JIT
    if (graph.empty()) return {};
    const int root = (rootNode >= 0) ? rootNode : graph.root();
    const uint64_t mkey = gpu_mem_key(graph, root);
    auto& mc = gpu_mem_cache();
    {
        std::lock_guard<std::mutex> lk(mc.mu);
        if (auto it = mc.ptx.find(mkey); it != mc.ptx.end()) {
            return it->second;
        }
    }

    const brass::cache::CacheKey dkey = build_gpu_artifact_key("bromesh.sdf.ptx/1", graph, root);
    auto& store = bromesh_artifact_cache();
    if (auto art = store.get(dkey)) {
        const std::string_view sv = art->text("ptx");
        if (!sv.empty()) {
            std::string out(sv);
            std::lock_guard<std::mutex> lk(mc.mu);
            mc.ptx[mkey] = out;
            return out;
        }
    }

    brass::Module mod("sdf_ptx_" + std::to_string(mkey));
    brass::Function* fn = buildGpuSdfKernel(mod, graph, root);
    std::string ptx = brass::target::PtxTarget::emit_function(*fn);
    if (!ptx.empty()) {
        std::string putErr;
        store.put(dkey, {brass::cache::BlobRef("ptx", ptx)}, &putErr);
        std::lock_guard<std::mutex> lk(mc.mu);
        mc.ptx[mkey] = ptx;
    }
    return ptx;
#else
    (void)graph;
    (void)rootNode;
    return {};
#endif
}

std::string JitSdfCompiler::compile_msl(const SdfGraph& graph, int rootNode) {
#if BROMESH_HAS_BRASS_JIT
    if (graph.empty()) return {};
    const int root = (rootNode >= 0) ? rootNode : graph.root();
    const uint64_t mkey = gpu_mem_key(graph, root);
    auto& mc = gpu_mem_cache();
    {
        std::lock_guard<std::mutex> lk(mc.mu);
        if (auto it = mc.msl.find(mkey); it != mc.msl.end()) {
            return it->second;
        }
    }

    const brass::cache::CacheKey dkey = build_gpu_artifact_key("bromesh.sdf.msl/1", graph, root);
    auto& store = bromesh_artifact_cache();
    if (auto art = store.get(dkey)) {
        const std::string_view sv = art->text("msl");
        if (!sv.empty()) {
            std::string out(sv);
            std::lock_guard<std::mutex> lk(mc.mu);
            mc.msl[mkey] = out;
            return out;
        }
    }

    brass::Module mod("sdf_msl_" + std::to_string(mkey));
    brass::Function* fn = buildGpuSdfKernel(mod, graph, root);
    brass::DiagnosticReporter portDiag;
    if (!brass::verify_gpu_portability(*fn, {}, &portDiag)) {
        return {};
    }

    brass::target::MslOptions mslOpts;
    mslOpts.workgroup_x = 64;
    auto msl = brass::target::MslTarget::compile(*fn, mslOpts);
    if (!msl.source.empty()) {
        std::string putErr;
        store.put(dkey, {brass::cache::BlobRef("msl", msl.source)}, &putErr);
        std::lock_guard<std::mutex> lk(mc.mu);
        mc.msl[mkey] = msl.source;
    }
    return msl.source;
#else
    (void)graph;
    (void)rootNode;
    return {};
#endif
}

bool JitSdfCompiler::export_kernel_bundle(const std::string& bundlePath, std::string* err) {
#if BROMESH_HAS_BRASS_JIT
    bromesh_artifact_cache().flush();
    brass::cache::BundleStats stats;
    return bromesh_artifact_cache().export_bundle(
        std::filesystem::path(bundlePath), "bromesh.sdf.bundle/1", &stats, err);
#else
    (void)bundlePath;
    if (err) *err = "BROMESH_HAS_BRASS_JIT is disabled";
    return false;
#endif
}

bool JitSdfCompiler::import_kernel_bundle(const std::string& bundlePath, std::string* err) {
#if BROMESH_HAS_BRASS_JIT
    brass::cache::BundleStats stats;
    std::string meta;
    return bromesh_artifact_cache().import_bundle(
        std::filesystem::path(bundlePath), &meta, &stats, err);
#else
    (void)bundlePath;
    if (err) *err = "BROMESH_HAS_BRASS_JIT is disabled";
    return false;
#endif
}

} // namespace bromesh
