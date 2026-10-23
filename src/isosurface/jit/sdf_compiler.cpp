#include "sdf_compiler_internal.h"

#if BROMESH_HAS_BRASS_JIT
#include <brass/mir/scalar_opt.hpp>
#include <brass/mir/smt_guard_elim.hpp>
#include <brass/mir/superopt.hpp>
#endif

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace bromesh {

#if BROMESH_HAS_BRASS_JIT
using namespace brass;
using namespace brass::codegen;

Value* lowerNoise3d(SdfKernelBuilder& kb, Value* px, Value* py, Value* pz) {
    Value* flr_x = kb.floor_f32(px);
    Value* flr_y = kb.floor_f32(py);
    Value* flr_z = kb.floor_f32(pz);

    Value* i0 = kb.builder().build_fptosi_i32_f32(flr_x);
    Value* j0 = kb.builder().build_fptosi_i32_f32(flr_y);
    Value* k0 = kb.builder().build_fptosi_i32_f32(flr_z);

    Value* one_i = kb.const_i32(1);
    Value* i1 = kb.builder().build_add(i0, one_i);
    Value* j1 = kb.builder().build_add(j0, one_i);
    Value* k1 = kb.builder().build_add(k0, one_i);

    Value* fx = kb.sub(px, flr_x);
    Value* fy = kb.sub(py, flr_y);
    Value* fz = kb.sub(pz, flr_z);

    Value* c_3 = kb.const_f32(3.0f);
    Value* c_2 = kb.const_f32(2.0f);
    Value* u = kb.mul(kb.mul(fx, fx), kb.sub(c_3, kb.mul(fx, c_2)));
    Value* v = kb.mul(kb.mul(fy, fy), kb.sub(c_3, kb.mul(fy, c_2)));
    Value* w = kb.mul(kb.mul(fz, fz), kb.sub(c_3, kb.mul(fz, c_2)));

    auto hash_corner = [&](Value* ix, Value* iy, Value* iz) -> Value* {
        Value* h = kb.builder().build_xor(ix, kb.builder().build_mul(iy, kb.const_i32(374761393)));
        h = kb.builder().build_xor(h, kb.builder().build_mul(iz, kb.const_i32(668265263)));
        Value* shift13 = kb.builder().build_lshr(h, kb.const_i32(13));
        h = kb.builder().build_mul(kb.builder().build_xor(h, shift13), kb.const_i32(1274126177));
        Value* shift16 = kb.builder().build_lshr(h, kb.const_i32(16));
        h = kb.builder().build_xor(h, shift16);

        Value* masked = kb.builder().build_and(h, kb.const_i32(0x7FFF));
        Value* h_f = kb.builder().build_sitofp_f32_i32(masked);

        Value* scale = kb.const_f32(2.0f / 32767.0f);
        Value* one_f = kb.const_f32(1.0f);
        return kb.sub(kb.mul(h_f, scale), one_f);
    };

    Value* h000 = hash_corner(i0, j0, k0);
    Value* h100 = hash_corner(i1, j0, k0);
    Value* h010 = hash_corner(i0, j1, k0);
    Value* h110 = hash_corner(i1, j1, k0);
    Value* h001 = hash_corner(i0, j0, k1);
    Value* h101 = hash_corner(i1, j0, k1);
    Value* h011 = hash_corner(i0, j1, k1);
    Value* h111 = hash_corner(i1, j1, k1);

    Value* x00 = kb.lerp_f32(h000, h100, u);
    Value* x10 = kb.lerp_f32(h010, h110, u);
    Value* x01 = kb.lerp_f32(h001, h101, u);
    Value* x11 = kb.lerp_f32(h011, h111, u);

    Value* y0 = kb.lerp_f32(x00, x10, v);
    Value* y1 = kb.lerp_f32(x01, x11, v);

    return kb.lerp_f32(y0, y1, w);
}

Value* lowerToMir(SdfKernelBuilder& kb, const SdfGraph& graph, int nodeId,
                  Value* px, Value* py, Value* pz, bool gpuTarget) {
    if (nodeId < 0 || nodeId >= static_cast<int>(graph.size())) {
        return kb.const_f32(1e9f);
    }
    const SdfNode& node = graph.node(nodeId);
    switch (node.op) {
        case SdfOp::Sphere: {
            Value* r = kb.const_f32(node.f0);
            return kb.sdf_sphere(px, py, pz, r);
        }
        case SdfOp::Box: {
            Value* bx = kb.const_f32(node.v0.x);
            Value* by = kb.const_f32(node.v0.y);
            Value* bz = kb.const_f32(node.v0.z);
            return kb.sdf_box(px, py, pz, bx, by, bz);
        }
        case SdfOp::RoundedBox: {
            Value* bx = kb.const_f32(node.v0.x);
            Value* by = kb.const_f32(node.v0.y);
            Value* bz = kb.const_f32(node.v0.z);
            Value* r = kb.const_f32(node.f0);
            return kb.sdf_rounded_box(px, py, pz, bx, by, bz, r);
        }
        case SdfOp::Cylinder: {
            Value* r = kb.const_f32(node.f0);
            Value* h = kb.const_f32(node.f1);
            return kb.sdf_cylinder(px, py, pz, r, h);
        }
        case SdfOp::Capsule: {
            Value* ax = kb.const_f32(node.v0.x);
            Value* ay = kb.const_f32(node.v0.y);
            Value* az = kb.const_f32(node.v0.z);
            Value* bx = kb.const_f32(node.v1.x);
            Value* by = kb.const_f32(node.v1.y);
            Value* bz = kb.const_f32(node.v1.z);
            Value* r = kb.const_f32(node.f0);
            return kb.sdf_capsule(px, py, pz, ax, ay, az, bx, by, bz, r);
        }
        case SdfOp::Torus: {
            Value* majR = kb.const_f32(node.f0);
            Value* minR = kb.const_f32(node.f1);
            return kb.sdf_torus(px, py, pz, majR, minR);
        }
        case SdfOp::Plane: {
            Value* nx = kb.const_f32(node.v0.x);
            Value* ny = kb.const_f32(node.v0.y);
            Value* nz = kb.const_f32(node.v0.z);
            Value* d = kb.const_f32(node.f0);
            return kb.sdf_plane(px, py, pz, nx, ny, nz, d);
        }
        case SdfOp::Union: {
            Value* d1 = lowerToMir(kb, graph, node.left, px, py, pz, gpuTarget);
            Value* d2 = lowerToMir(kb, graph, node.right, px, py, pz, gpuTarget);
            return kb.op_union(d1, d2);
        }
        case SdfOp::Intersection: {
            Value* d1 = lowerToMir(kb, graph, node.left, px, py, pz, gpuTarget);
            Value* d2 = lowerToMir(kb, graph, node.right, px, py, pz, gpuTarget);
            return kb.op_intersection(d1, d2);
        }
        case SdfOp::Subtraction: {
            Value* d1 = lowerToMir(kb, graph, node.left, px, py, pz, gpuTarget);
            Value* d2 = lowerToMir(kb, graph, node.right, px, py, pz, gpuTarget);
            return kb.op_subtraction(d1, d2);
        }
        case SdfOp::SmoothUnion: {
            Value* d1 = lowerToMir(kb, graph, node.left, px, py, pz, gpuTarget);
            Value* d2 = lowerToMir(kb, graph, node.right, px, py, pz, gpuTarget);
            Value* k = kb.const_f32(node.f0);
            return kb.op_smooth_union(d1, d2, k);
        }
        case SdfOp::SmoothIntersection: {
            Value* d1 = lowerToMir(kb, graph, node.left, px, py, pz, gpuTarget);
            Value* d2 = lowerToMir(kb, graph, node.right, px, py, pz, gpuTarget);
            Value* k = kb.const_f32(node.f0);
            return kb.op_smooth_intersection(d1, d2, k);
        }
        case SdfOp::SmoothSubtraction: {
            Value* d1 = lowerToMir(kb, graph, node.left, px, py, pz, gpuTarget);
            Value* d2 = lowerToMir(kb, graph, node.right, px, py, pz, gpuTarget);
            Value* k = kb.const_f32(node.f0);
            return kb.op_smooth_subtraction(d1, d2, k);
        }
        case SdfOp::Translate: {
            Value* ox = kb.const_f32(node.v0.x);
            Value* oy = kb.const_f32(node.v0.y);
            Value* oz = kb.const_f32(node.v0.z);
            Value* tx = nullptr, *ty = nullptr, *tz = nullptr;
            kb.translate(px, py, pz, ox, oy, oz, tx, ty, tz);
            return lowerToMir(kb, graph, node.left, tx, ty, tz, gpuTarget);
        }
        case SdfOp::RotateY: {
            Value* c = kb.const_f32(std::cos(node.f0));
            Value* s = kb.const_f32(std::sin(node.f0));
            Value* rx = kb.sub(kb.mul(px, c), kb.mul(pz, s));
            Value* rz = kb.add(kb.mul(px, s), kb.mul(pz, c));
            return lowerToMir(kb, graph, node.left, rx, py, rz, gpuTarget);
        }
        case SdfOp::ScaleUniform: {
            Value* s = kb.const_f32(node.f0);
            Value* sx = nullptr, *sy = nullptr, *sz = nullptr;
            kb.scale_uniform(px, py, pz, s, sx, sy, sz);
            Value* dist = lowerToMir(kb, graph, node.left, sx, sy, sz, gpuTarget);
            return kb.scale_uniform_dist(dist, s);
        }
        case SdfOp::Noise3D: {
            Value* nx = kb.mul(px, kb.const_f32(node.v0.x));
            Value* ny = kb.mul(py, kb.const_f32(node.v0.y));
            Value* nz = kb.mul(pz, kb.const_f32(node.v0.z));
            Value* n = lowerNoise3d(kb, nx, ny, nz);
            return kb.mul(n, kb.const_f32(node.f0));
        }
        case SdfOp::Displace: {
            Value* d = lowerToMir(kb, graph, node.left, px, py, pz, gpuTarget);
            Value* n = lowerToMir(kb, graph, node.right, px, py, pz, gpuTarget);
            return kb.add(d, n);
        }
        case SdfOp::CustomExpr:
            return lowerCustomExprToMir(kb, graph.customExprData(node.left), px, py, pz, gpuTarget);
    }
    return kb.const_f32(1e9f);
}

static Function* buildGpuSdfKernel(Module& mod, const SdfGraph& graph, int rootNode) {
    int root = rootNode >= 0 ? rootNode : graph.root();
    Function* fn = mod.create_function("eval_sdf_gpu", Type::void_type(), {
        Type::ptr(), // float* out_grid
        Type::i32(), // int32_t dimX
        Type::i32(), // int32_t dimY
        Type::i32(), // int32_t dimZ
        Type::f32(), // float minX
        Type::f32(), // float minY
        Type::f32(), // float minZ
        Type::f32(), // float stepX
        Type::f32(), // float stepY
        Type::f32()  // float stepZ
    });

    SdfKernelBuilder kb(mod, fn);
    BasicBlock* entry = kb.builder().append_block("entry");
    Value* out_grid = kb.builder().add_block_param(entry, Type::ptr());
    Value* dimX = kb.builder().add_block_param(entry, Type::i32());
    Value* dimY = kb.builder().add_block_param(entry, Type::i32());
    Value* dimZ = kb.builder().add_block_param(entry, Type::i32());
    Value* minX = kb.builder().add_block_param(entry, Type::f32());
    Value* minY = kb.builder().add_block_param(entry, Type::f32());
    Value* minZ = kb.builder().add_block_param(entry, Type::f32());
    Value* stepX = kb.builder().add_block_param(entry, Type::f32());
    Value* stepY = kb.builder().add_block_param(entry, Type::f32());
    Value* stepZ = kb.builder().add_block_param(entry, Type::f32());

    kb.position_at_end(entry);
    Value* idx = kb.global_tid_x();
    Value* stride_z = kb.mul(dimX, dimY);
    Value* total = kb.mul(stride_z, dimZ);
    Value* in_bounds = kb.builder().build_slt(idx, total);

    kb.if_then(in_bounds, [&]() {
        Value* iz = kb.div(idx, stride_z);
        Value* rem = kb.sub(idx, kb.mul(iz, stride_z));
        Value* iy = kb.div(rem, dimX);
        Value* ix = kb.sub(rem, kb.mul(iy, dimX));

        Value* px = kb.fma(kb.to_f32(ix), stepX, minX);
        Value* py = kb.fma(kb.to_f32(iy), stepY, minY);
        Value* pz = kb.fma(kb.to_f32(iz), stepZ, minZ);

        Value* sdf_val = lowerToMir(kb, graph, root, px, py, pz, true);
        kb.store_f32_indexed(out_grid, idx, sdf_val, 4, 0);
    });
    kb.builder().build_ret_void();

    run_smt_guard_elimination(*fn, mod);
    run_mir_superopt(*fn, mod);
    eliminate_dead_code(*fn);

    return fn;
}
#endif

float evalScalarNode(const SdfGraph& graph, int nodeId, float px, float py, float pz) {
    if (nodeId < 0 || nodeId >= static_cast<int>(graph.size())) {
        return 1e9f;
    }
    const SdfNode& node = graph.node(nodeId);
    switch (node.op) {
        case SdfOp::Sphere:
            return std::sqrt(px * px + py * py + pz * pz) - node.f0;
        case SdfOp::Box: {
            float qx = std::fabs(px) - node.v0.x;
            float qy = std::fabs(py) - node.v0.y;
            float qz = std::fabs(pz) - node.v0.z;
            float ox = std::max(qx, 0.0f), oy = std::max(qy, 0.0f), oz = std::max(qz, 0.0f);
            return std::sqrt(ox * ox + oy * oy + oz * oz) + std::min(std::max({qx, qy, qz}), 0.0f);
        }
        case SdfOp::RoundedBox: {
            float r = node.f0;
            float qx = std::fabs(px) - node.v0.x + r;
            float qy = std::fabs(py) - node.v0.y + r;
            float qz = std::fabs(pz) - node.v0.z + r;
            float ox = std::max(qx, 0.0f), oy = std::max(qy, 0.0f), oz = std::max(qz, 0.0f);
            return std::sqrt(ox * ox + oy * oy + oz * oz) + std::min(std::max({qx, qy, qz}), 0.0f) - r;
        }
        case SdfOp::Cylinder: {
            float dx = std::sqrt(px * px + pz * pz) - node.f0;
            float dy = std::fabs(py) - node.f1;
            float ox = std::max(dx, 0.0f), oy = std::max(dy, 0.0f);
            return std::sqrt(ox * ox + oy * oy) + std::min(std::max(dx, dy), 0.0f);
        }
        case SdfOp::Capsule: {
            float pax = px - node.v0.x, pay = py - node.v0.y, paz = pz - node.v0.z;
            float bax = node.v1.x - node.v0.x, bay = node.v1.y - node.v0.y, baz = node.v1.z - node.v0.z;
            float ba_len2 = bax * bax + bay * bay + baz * baz;
            float h = (ba_len2 > 1e-12f) ? std::clamp((pax * bax + pay * bay + paz * baz) / ba_len2, 0.0f, 1.0f) : 0.0f;
            float rx = pax - bax * h, ry = pay - bay * h, rz = paz - baz * h;
            return std::sqrt(rx * rx + ry * ry + rz * rz) - node.f0;
        }
        case SdfOp::Torus: {
            float qx = std::sqrt(px * px + pz * pz) - node.f0;
            return std::sqrt(qx * qx + py * py) - node.f1;
        }
        case SdfOp::Plane:
            return (px * node.v0.x + py * node.v0.y + pz * node.v0.z) + node.f0;
        case SdfOp::Union:
            return std::min(evalScalarNode(graph, node.left, px, py, pz),
                            evalScalarNode(graph, node.right, px, py, pz));
        case SdfOp::Intersection:
            return std::max(evalScalarNode(graph, node.left, px, py, pz),
                            evalScalarNode(graph, node.right, px, py, pz));
        case SdfOp::Subtraction:
            return std::max(evalScalarNode(graph, node.left, px, py, pz),
                            -evalScalarNode(graph, node.right, px, py, pz));
        case SdfOp::SmoothUnion: {
            float d1 = evalScalarNode(graph, node.left, px, py, pz);
            float d2 = evalScalarNode(graph, node.right, px, py, pz);
            float h = std::clamp(0.5f + 0.5f * (d2 - d1) / node.f0, 0.0f, 1.0f);
            return (d2 + (d1 - d2) * h) - node.f0 * h * (1.0f - h);
        }
        case SdfOp::SmoothIntersection: {
            float d1 = evalScalarNode(graph, node.left, px, py, pz);
            float d2 = evalScalarNode(graph, node.right, px, py, pz);
            float h = std::clamp(0.5f - 0.5f * (d2 - d1) / node.f0, 0.0f, 1.0f);
            return (d2 + (d1 - d2) * h) + node.f0 * h * (1.0f - h);
        }
        case SdfOp::SmoothSubtraction: {
            float d1 = evalScalarNode(graph, node.left, px, py, pz);
            float d2 = evalScalarNode(graph, node.right, px, py, pz);
            float h = std::clamp(0.5f - 0.5f * (d2 + d1) / node.f0, 0.0f, 1.0f);
            return (d1 + (-d2 - d1) * h) + node.f0 * h * (1.0f - h);
        }
        case SdfOp::Translate:
            return evalScalarNode(graph, node.left, px - node.v0.x, py - node.v0.y, pz - node.v0.z);
        case SdfOp::RotateY: {
            float c = std::cos(node.f0), s = std::sin(node.f0);
            return evalScalarNode(graph, node.left, px * c - pz * s, py, px * s + pz * c);
        }
        case SdfOp::ScaleUniform: {
            float s = node.f0;
            return evalScalarNode(graph, node.left, px / s, py / s, pz / s) * s;
        }
        case SdfOp::Noise3D: {
            float nx = px * node.v0.x, ny = py * node.v0.y, nz = pz * node.v0.z;
            float flrx = std::floor(nx), flry = std::floor(ny), flrz = std::floor(nz);
            int32_t i0 = static_cast<int32_t>(flrx), j0 = static_cast<int32_t>(flry), k0 = static_cast<int32_t>(flrz);
            float fx = nx - flrx, fy = ny - flry, fz = nz - flrz;
            float u = fx * fx * (3.0f - 2.0f * fx), v = fy * fy * (3.0f - 2.0f * fy), w = fz * fz * (3.0f - 2.0f * fz);
            auto hashc = [](int32_t ix, int32_t iy, int32_t iz) -> float {
                uint32_t h = static_cast<uint32_t>(ix) ^ (static_cast<uint32_t>(iy) * 374761393u);
                h ^= static_cast<uint32_t>(iz) * 668265263u;
                h = (h ^ (h >> 13)) * 1274126177u;
                h ^= h >> 16;
                return static_cast<float>(h & 0x7FFFu) * (2.0f / 32767.0f) - 1.0f;
            };
            float h000 = hashc(i0, j0, k0), h100 = hashc(i0 + 1, j0, k0);
            float h010 = hashc(i0, j0 + 1, k0), h110 = hashc(i0 + 1, j0 + 1, k0);
            float h001 = hashc(i0, j0, k0 + 1), h101 = hashc(i0 + 1, j0, k0 + 1);
            float h011 = hashc(i0, j0 + 1, k0 + 1), h111 = hashc(i0 + 1, j0 + 1, k0 + 1);
            float x00 = h000 + (h100 - h000) * u, x10 = h010 + (h110 - h010) * u;
            float x01 = h001 + (h101 - h001) * u, x11 = h011 + (h111 - h011) * u;
            float y0 = x00 + (x10 - x00) * v, y1 = x01 + (x11 - x01) * v;
            return (y0 + (y1 - y0) * w) * node.f0;
        }
        case SdfOp::Displace:
            return evalScalarNode(graph, node.left, px, py, pz) + evalScalarNode(graph, node.right, px, py, pz);
        case SdfOp::CustomExpr:
            return evalCustomExprScalar(graph.customExprData(node.left), px, py, pz);
    }
    return 1e9f;
}

JitSdfCompiler& JitSdfCompiler::instance() {
    static JitSdfCompiler s_instance;
    return s_instance;
}

void JitSdfCompiler::clearCache() {
    std::lock_guard<std::mutex> lock(m_cacheMutex);
    m_cache.clear();
}

CompiledSdfKernel JitSdfCompiler::compile(const SdfGraph& graph, int rootNode) {
    int root = rootNode >= 0 ? rootNode : graph.root();
    uint64_t hash = graph.computeHash();
    if (rootNode >= 0 && rootNode != graph.root()) {
        hash ^= (static_cast<uint64_t>(rootNode) * 0x9e3779b97f4a7c15ULL);
    }

    {
        std::lock_guard<std::mutex> lock(m_cacheMutex);
        auto it = m_cache.find(hash);
        if (it != m_cache.end()) {
            return it->second;
        }
    }

    CompiledSdfKernel kernel;
    kernel.hash = hash;

#if BROMESH_HAS_BRASS_JIT
    if (Target::host().is_x64() || Target::host().is_aarch64()) {
        try {
            KernelOptions opts;
            opts.opt_level = 2;
            opts.enable_optimizations = true;
            opts.enable_sroa = true;
            opts.enable_bce = true;
            opts.enable_smt_guard_elim = true;
            opts.enable_superopt = true;
            opts.enable_avx2 = true;
            opts.enable_fma = true;
            opts.enable_parallel = false;

            // 1. Grid Evaluation Function
            Module modGrid("sdf_grid_" + std::to_string(hash));
            Function* fnGrid = modGrid.create_function("eval_grid", Type::void_type(), {
                Type::ptr(), Type::ptr(), Type::ptr(), Type::i32(), Type::i32(), Type::i32()
            });

            SdfKernelBuilder kbGrid(modGrid, fnGrid);
            BasicBlock* entryGrid = kbGrid.builder().append_block("entry");
            Value* outField = kbGrid.builder().add_block_param(entryGrid, Type::ptr());
            Value* boundsMin = kbGrid.builder().add_block_param(entryGrid, Type::ptr());
            Value* cellSize = kbGrid.builder().add_block_param(entryGrid, Type::ptr());
            Value* dimX = kbGrid.builder().add_block_param(entryGrid, Type::i32());
            Value* dimY = kbGrid.builder().add_block_param(entryGrid, Type::i32());
            Value* dimZ = kbGrid.builder().add_block_param(entryGrid, Type::i32());

            kbGrid.position_at_end(entryGrid);
            Value* dimX64 = kbGrid.builder().build_sext_i64(dimX);
            Value* dimY64 = kbGrid.builder().build_sext_i64(dimY);
            Value* dimZ64 = kbGrid.builder().build_sext_i64(dimZ);

            Value* ox = kbGrid.load_f32(boundsMin, 0);
            Value* oy = kbGrid.load_f32(boundsMin, 4);
            Value* oz = kbGrid.load_f32(boundsMin, 8);
            Value* sx = kbGrid.load_f32(cellSize, 0);
            Value* sy = kbGrid.load_f32(cellSize, 4);
            Value* sz = kbGrid.load_f32(cellSize, 8);

            kbGrid.emit_grid_eval_loop(outField, dimX64, dimY64, dimZ64, ox, oy, oz, sx, sy, sz,
                                      [&](Value* px, Value* py, Value* pz) {
                return lowerToMir(kbGrid, graph, root, px, py, pz, false);
            });
            kbGrid.builder().build_ret_void();

            KernelJit jitGrid(opts);
            kernel.gridKernel = jitGrid.compile(*fnGrid);
            if (kernel.gridKernel.is_valid()) {
                kernel.gridFn = kernel.gridKernel.as<SdfGridEvalFn>();
            }

            // 1b. Span Evaluation Function (bit-exact global X indexing for macro-block pruning)
            Module modSpan("sdf_span_" + std::to_string(hash));
            Function* fnSpan = modSpan.create_function("eval_span", Type::void_type(), {
                Type::ptr(), Type::f32(), Type::f32(), Type::f32(), Type::f32(), Type::i64(), Type::i64()
            });
            SdfKernelBuilder kbSpan(modSpan, fnSpan);
            BasicBlock* entrySpan = kbSpan.builder().append_block("entry");
            Value* outRow = kbSpan.builder().add_block_param(entrySpan, Type::ptr());
            Value* oxSpan = kbSpan.builder().add_block_param(entrySpan, Type::f32());
            Value* pySpan = kbSpan.builder().add_block_param(entrySpan, Type::f32());
            Value* pzSpan = kbSpan.builder().add_block_param(entrySpan, Type::f32());
            Value* sxSpan = kbSpan.builder().add_block_param(entrySpan, Type::f32());
            Value* xStart = kbSpan.builder().add_block_param(entrySpan, Type::i64());
            Value* xEnd = kbSpan.builder().add_block_param(entrySpan, Type::i64());

            kbSpan.position_at_end(entrySpan);
            kbSpan.for_range(xStart, xEnd, kbSpan.const_i64(1), [&](Value* x) {
                Value* px = kbSpan.fma(kbSpan.to_f32(x), sxSpan, oxSpan);
                Value* dist = lowerToMir(kbSpan, graph, root, px, pySpan, pzSpan, false);
                kbSpan.store_f32_indexed(outRow, x, dist, 4, 0);
            });
            kbSpan.builder().build_ret_void();

            KernelJit jitSpan(opts);
            kernel.spanKernel = jitSpan.compile(*fnSpan);
            if (kernel.spanKernel.is_valid()) {
                kernel.spanFn = kernel.spanKernel.as<SdfSpanEvalFn>();
            }

            // 2. Point Evaluation Function
            Module modPoint("sdf_point_" + std::to_string(hash));
            Function* fnPoint = modPoint.create_function("eval_point", Type::f32(), {
                Type::f32(), Type::f32(), Type::f32()
            });

            SdfKernelBuilder kbPoint(modPoint, fnPoint);
            BasicBlock* entryPoint = kbPoint.builder().append_block("entry");
            Value* px = kbPoint.builder().add_block_param(entryPoint, Type::f32());
            Value* py = kbPoint.builder().add_block_param(entryPoint, Type::f32());
            Value* pz = kbPoint.builder().add_block_param(entryPoint, Type::f32());

            kbPoint.position_at_end(entryPoint);
            Value* d = lowerToMir(kbPoint, graph, root, px, py, pz, false);
            kbPoint.builder().build_ret(d);

            KernelJit jitPoint(opts);
            kernel.pointKernel = jitPoint.compile(*fnPoint);
            if (kernel.pointKernel.is_valid()) {
                kernel.pointFn = kernel.pointKernel.as<SdfPointEvalFn>();
            }

            // 3. 1-Pass Symbolic Auto-Diff Gradient Function (#107)
            Module modGrad("sdf_grad_" + std::to_string(hash));
            Function* fnGrad = buildSdfGradFunction(modGrad, graph, root);
            KernelJit jitGrad(opts);
            kernel.gradKernel = jitGrad.compile(*fnGrad);
            if (kernel.gradKernel.is_valid()) {
                kernel.gradFn = kernel.gradKernel.as<SdfGradEvalFn>();
            }

            // 4. Interval / Lipschitz Macro-Block Bounds Function (#106)
            Module modBounds("sdf_bounds_" + std::to_string(hash));
            Function* fnBounds = buildSdfBlockBoundsFunction(modBounds, graph, root);
            KernelJit jitBounds(opts);
            kernel.blockBoundsKernel = jitBounds.compile(*fnBounds);
            if (kernel.blockBoundsKernel.is_valid()) {
                kernel.blockBoundsFn = kernel.blockBoundsKernel.as<SdfBlockBoundsFn>();
            }

            kernel.isJit = (kernel.gridFn != nullptr && kernel.pointFn != nullptr &&
                            kernel.gradFn != nullptr && kernel.blockBoundsFn != nullptr);
        } catch (const std::exception& e) {
            std::cerr << "JitSdfCompiler compilation failed: " << e.what() << std::endl;
        }
    }
#endif

    {
        std::lock_guard<std::mutex> lock(m_cacheMutex);
        m_cache[hash] = kernel;
    }
    return kernel;
}

std::vector<uint32_t> JitSdfCompiler::compile_spirv(const SdfGraph& graph, int rootNode) {
#if BROMESH_HAS_BRASS_JIT
    uint64_t hash = graph.computeHash();
    Module mod("sdf_spirv_" + std::to_string(hash));
    Function* fn = buildGpuSdfKernel(mod, graph, rootNode);
    auto kernel = brass::target::SpirvTarget::compile(*fn);
    return kernel.words;
#else
    (void)graph;
    (void)rootNode;
    return {};
#endif
}

std::string JitSdfCompiler::compile_ptx(const SdfGraph& graph, int rootNode) {
#if BROMESH_HAS_BRASS_JIT
    uint64_t hash = graph.computeHash();
    Module mod("sdf_ptx_" + std::to_string(hash));
    Function* fn = buildGpuSdfKernel(mod, graph, rootNode);
    return brass::target::PtxTarget::emit_function(*fn);
#else
    (void)graph;
    (void)rootNode;
    return {};
#endif
}

SDFVolume JitSdfCompiler::evaluateVolume(const SdfGraph& graph, int dimX, int dimY, int dimZ,
                                        const bromath::AABB3& bounds, int rootNode, bool parallel,
                                        float isoLevel, bool enableMacroBlockPruning) {
    SDFVolume vol;
    if (dimX <= 0 || dimY <= 0 || dimZ <= 0) {
        return vol;
    }
    vol.dimX = dimX;
    vol.dimY = dimY;
    vol.dimZ = dimZ;
    vol.bounds = bounds;

    vol.cellSize[0] = dimX > 1 ? (bounds.max.x - bounds.min.x) / static_cast<float>(dimX - 1) : 1.0f;
    vol.cellSize[1] = dimY > 1 ? (bounds.max.y - bounds.min.y) / static_cast<float>(dimY - 1) : 1.0f;
    vol.cellSize[2] = dimZ > 1 ? (bounds.max.z - bounds.min.z) / static_cast<float>(dimZ - 1) : 1.0f;

    const size_t totalSamples = static_cast<size_t>(dimX) * dimY * dimZ;
    vol.field.resize(totalSamples, 0.0f);

    int root = rootNode >= 0 ? rootNode : graph.root();
    CompiledSdfKernel kernel = compile(graph, root);

    const int cellsX = std::max(1, dimX - 1);
    const int cellsY = std::max(1, dimY - 1);
    const int cellsZ = std::max(1, dimZ - 1);
    const int blocksX = (cellsX + 7) / 8;
    const int blocksY = (cellsY + 7) / 8;
    const int blocksZ = (cellsZ + 7) / 8;
    const size_t totalBlocks = static_cast<size_t>(blocksX) * blocksY * blocksZ;

    vol.blocksX = blocksX;
    vol.blocksY = blocksY;
    vol.blocksZ = blocksZ;
    vol.totalBlockCount = totalBlocks;

    std::vector<float> blockFillVal;
    bool usePrunedEval = false;

    if (enableMacroBlockPruning && blocksX >= 2 && blocksY >= 2 && blocksZ >= 2) {
        vol.activeBlockMask.assign(totalBlocks, 0);
        blockFillVal.resize(totalBlocks, 0.0f);

        auto evalBlockSliceRange = [&](int64_t bzStart, int64_t bzEnd) {
            const float margin = 1e-4f;
            for (int64_t bz = bzStart; bz < bzEnd; ++bz) {
                int gz0 = static_cast<int>(bz) * 8;
                int gz1 = std::min(gz0 + 8, cellsZ);
                float minZ = bounds.min.z + static_cast<float>(gz0) * vol.cellSize[2];
                float maxZ = bounds.min.z + static_cast<float>(gz1) * vol.cellSize[2];

                for (int by = 0; by < blocksY; ++by) {
                    int gy0 = by * 8;
                    int gy1 = std::min(gy0 + 8, cellsY);
                    float minY = bounds.min.y + static_cast<float>(gy0) * vol.cellSize[1];
                    float maxY = bounds.min.y + static_cast<float>(gy1) * vol.cellSize[1];

                    for (int bx = 0; bx < blocksX; ++bx) {
                        int gx0 = bx * 8;
                        int gx1 = std::min(gx0 + 8, cellsX);
                        float minX = bounds.min.x + static_cast<float>(gx0) * vol.cellSize[0];
                        float maxX = bounds.min.x + static_cast<float>(gx1) * vol.cellSize[0];

                        float dMin = 0.0f, dMax = 0.0f;
                        if (kernel.blockBoundsFn) {
                            kernel.blockBoundsFn(minX, minY, minZ, maxX, maxY, maxZ, &dMin, &dMax);
                        } else {
                            ::bromesh::evalBlockBounds(graph, root, {{minX, minY, minZ}, {maxX, maxY, maxZ}}, dMin, dMax);
                        }

                        bool isBoundary = (bx == 0 || bx == blocksX - 1 ||
                                           by == 0 || by == blocksY - 1 ||
                                           bz == 0 || bz == blocksZ - 1);
                        size_t bIdx = (static_cast<size_t>(bz) * blocksY + by) * blocksX + bx;
                        if (dMin > isoLevel + margin) {
                            vol.activeBlockMask[bIdx] = 0;
                            blockFillVal[bIdx] = dMin;
                        } else if (dMax < isoLevel - margin && !isBoundary) {
                            vol.activeBlockMask[bIdx] = 0;
                            blockFillVal[bIdx] = dMax;
                        } else {
                            vol.activeBlockMask[bIdx] = 1;
                            blockFillVal[bIdx] = 0.0f;
                        }
                    }
                }
            }
        };

#if BROMESH_HAS_BRASS_JIT
        if (parallel && blocksZ >= 4) {
            brass::codegen::parallel_for(0, static_cast<int64_t>(blocksZ), evalBlockSliceRange);
        } else {
            evalBlockSliceRange(0, static_cast<int64_t>(blocksZ));
        }
#else
        evalBlockSliceRange(0, static_cast<int64_t>(blocksZ));
#endif

        size_t activeCount = 0;
        for (uint8_t m : vol.activeBlockMask) activeCount += m;
        vol.activeBlockCount = activeCount;
        usePrunedEval = (activeCount < totalBlocks);
    } else {
        vol.activeBlockCount = totalBlocks;
    }

    if (kernel.gridFn) {
        auto evalSliceRange = [&](int64_t zStart, int64_t zEnd) {
            if (!usePrunedEval) {
                for (int64_t z = zStart; z < zEnd; ++z) {
                    float sliceBoundsMin[3] = {
                        bounds.min.x,
                        bounds.min.y,
                        bounds.min.z + static_cast<float>(z) * vol.cellSize[2]
                    };
                    float* sliceOut = vol.field.data() + static_cast<size_t>(z) * static_cast<size_t>(dimY) * static_cast<size_t>(dimX);
                    kernel.gridFn(sliceOut, sliceBoundsMin, vol.cellSize, dimX, dimY, 1);
                }
                return;
            }

            for (int64_t z = zStart; z < zEnd; ++z) {
                int bz0 = (z > 0) ? static_cast<int>((z - 1) >> 3) : 0;
                int bz1 = std::min(static_cast<int>(z >> 3), blocksZ - 1);
                float pz = bounds.min.z + static_cast<float>(z) * vol.cellSize[2];

                for (int y = 0; y < dimY; ++y) {
                    int by0 = (y > 0) ? ((y - 1) >> 3) : 0;
                    int by1 = std::min(y >> 3, blocksY - 1);
                    float py = std::fma(static_cast<float>(y), vol.cellSize[1], bounds.min.y);
                    float* rowOut = vol.field.data() + (static_cast<size_t>(z) * dimY + y) * dimX;

                    auto isColActive = [&](int bx) -> bool {
                        return vol.activeBlockMask[(static_cast<size_t>(bz0) * blocksY + by0) * blocksX + bx] != 0 ||
                               vol.activeBlockMask[(static_cast<size_t>(bz0) * blocksY + by1) * blocksX + bx] != 0 ||
                               vol.activeBlockMask[(static_cast<size_t>(bz1) * blocksY + by0) * blocksX + bx] != 0 ||
                               vol.activeBlockMask[(static_cast<size_t>(bz1) * blocksY + by1) * blocksX + bx] != 0;
                    };

                    int nextFillX = 0;
                    int bx = 0;
                    while (bx < blocksX) {
                        if (!isColActive(bx)) {
                            int fillEnd = std::min((bx + 1) * 8, dimX);
                            if (bx + 1 < blocksX && isColActive(bx + 1)) {
                                fillEnd = (bx + 1) * 8;
                            }
                            if (fillEnd > nextFillX) {
                                float fillVal = blockFillVal[(static_cast<size_t>(bz1) * blocksY + by1) * blocksX + bx];
                                std::fill(rowOut + nextFillX, rowOut + fillEnd, fillVal);
                                nextFillX = fillEnd;
                            }
                            ++bx;
                        } else {
                            int bxStart = bx;
                            while (bx + 1 < blocksX && isColActive(bx + 1)) {
                                ++bx;
                            }
                            int gxStart = bxStart * 8;
                            int gxEnd = std::min((bx + 1) * 8, dimX - 1);
                            if (kernel.spanFn) {
                                kernel.spanFn(rowOut, bounds.min.x, py, pz, vol.cellSize[0],
                                              static_cast<int64_t>(gxStart), static_cast<int64_t>(gxEnd + 1));
                            } else {
                                int gxCount = gxEnd - gxStart + 1;
                                float segBoundsMin[3] = {
                                    bounds.min.x + static_cast<float>(gxStart) * vol.cellSize[0],
                                    py,
                                    pz
                                };
                                kernel.gridFn(rowOut + gxStart, segBoundsMin, vol.cellSize, gxCount, 1, 1);
                            }
                            nextFillX = gxEnd + 1;
                            ++bx;
                        }
                    }
                }
            }
        };
#if BROMESH_HAS_BRASS_JIT
        if (parallel && dimZ >= 2) {
            brass::codegen::parallel_for(0, static_cast<int64_t>(dimZ), evalSliceRange);
        } else {
            evalSliceRange(0, static_cast<int64_t>(dimZ));
        }
#else
        (void)parallel;
        evalSliceRange(0, static_cast<int64_t>(dimZ));
#endif
    } else {
        #pragma omp parallel for
        for (int z = 0; z < dimZ; ++z) {
            float pz = bounds.min.z + z * vol.cellSize[2];
            for (int y = 0; y < dimY; ++y) {
                float py = bounds.min.y + y * vol.cellSize[1];
                for (int x = 0; x < dimX; ++x) {
                    float px = bounds.min.x + x * vol.cellSize[0];
                    size_t idx = (static_cast<size_t>(z) * dimY + y) * dimX + x;
                    vol.field[idx] = evalScalarNode(graph, root, px, py, pz);
                }
            }
        }
    }

    return vol;
}

} // namespace bromesh
