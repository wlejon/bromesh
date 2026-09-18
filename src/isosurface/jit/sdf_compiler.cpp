#include "bromesh/isosurface/jit/sdf_compiler.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace bromesh {

#if BROMESH_HAS_BRASS_JIT
using namespace brass;
using namespace brass::codegen;

static Value* lowerToMir(SdfKernelBuilder& kb, const SdfGraph& graph, int nodeId,
                         Value* px, Value* py, Value* pz) {
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
            Value* d1 = lowerToMir(kb, graph, node.left, px, py, pz);
            Value* d2 = lowerToMir(kb, graph, node.right, px, py, pz);
            return kb.op_union(d1, d2);
        }
        case SdfOp::Intersection: {
            Value* d1 = lowerToMir(kb, graph, node.left, px, py, pz);
            Value* d2 = lowerToMir(kb, graph, node.right, px, py, pz);
            return kb.op_intersection(d1, d2);
        }
        case SdfOp::Subtraction: {
            Value* d1 = lowerToMir(kb, graph, node.left, px, py, pz);
            Value* d2 = lowerToMir(kb, graph, node.right, px, py, pz);
            return kb.op_subtraction(d1, d2);
        }
        case SdfOp::SmoothUnion: {
            Value* d1 = lowerToMir(kb, graph, node.left, px, py, pz);
            Value* d2 = lowerToMir(kb, graph, node.right, px, py, pz);
            Value* k = kb.const_f32(node.f0);
            return kb.op_smooth_union(d1, d2, k);
        }
        case SdfOp::SmoothIntersection: {
            Value* d1 = lowerToMir(kb, graph, node.left, px, py, pz);
            Value* d2 = lowerToMir(kb, graph, node.right, px, py, pz);
            Value* k = kb.const_f32(node.f0);
            return kb.op_smooth_intersection(d1, d2, k);
        }
        case SdfOp::SmoothSubtraction: {
            Value* d1 = lowerToMir(kb, graph, node.left, px, py, pz);
            Value* d2 = lowerToMir(kb, graph, node.right, px, py, pz);
            Value* k = kb.const_f32(node.f0);
            return kb.op_smooth_subtraction(d1, d2, k);
        }
        case SdfOp::Translate: {
            Value* ox = kb.const_f32(node.v0.x);
            Value* oy = kb.const_f32(node.v0.y);
            Value* oz = kb.const_f32(node.v0.z);
            Value* tx = nullptr;
            Value* ty = nullptr;
            Value* tz = nullptr;
            kb.translate(px, py, pz, ox, oy, oz, tx, ty, tz);
            return lowerToMir(kb, graph, node.left, tx, ty, tz);
        }
        case SdfOp::RotateY: {
            Value* angle = kb.const_f32(node.f0);
            Value* rx = nullptr;
            Value* ry = nullptr;
            Value* rz = nullptr;
            kb.rotate_y(px, py, pz, angle, rx, ry, rz);
            return lowerToMir(kb, graph, node.left, rx, ry, rz);
        }
        case SdfOp::ScaleUniform: {
            Value* s = kb.const_f32(node.f0);
            Value* sx = nullptr;
            Value* sy = nullptr;
            Value* sz = nullptr;
            kb.scale_uniform(px, py, pz, s, sx, sy, sz);
            Value* dist = lowerToMir(kb, graph, node.left, sx, sy, sz);
            return kb.scale_uniform_dist(dist, s);
        }
        case SdfOp::Noise3D: {
            Value* nx = kb.mul(px, kb.const_f32(node.v0.x));
            Value* ny = kb.mul(py, kb.const_f32(node.v0.y));
            Value* nz = kb.mul(pz, kb.const_f32(node.v0.z));
            Value* n = kb.noise_3d(nx, ny, nz);
            return kb.mul(n, kb.const_f32(node.f0));
        }
        case SdfOp::Displace: {
            Value* d = lowerToMir(kb, graph, node.left, px, py, pz);
            Value* n = lowerToMir(kb, graph, node.right, px, py, pz);
            return kb.add(d, n);
        }
    }
    return kb.const_f32(1e9f);
}
#endif

static float evalScalarNode(const SdfGraph& graph, int nodeId, float px, float py, float pz) {
    if (nodeId < 0 || nodeId >= static_cast<int>(graph.size())) {
        return 1e9f;
    }
    const SdfNode& node = graph.node(nodeId);
    switch (node.op) {
        case SdfOp::Sphere: {
            return std::sqrt(px * px + py * py + pz * pz) - node.f0;
        }
        case SdfOp::Box: {
            float qx = std::fabs(px) - node.v0.x;
            float qy = std::fabs(py) - node.v0.y;
            float qz = std::fabs(pz) - node.v0.z;
            float ox = std::max(qx, 0.0f);
            float oy = std::max(qy, 0.0f);
            float oz = std::max(qz, 0.0f);
            float outside = std::sqrt(ox * ox + oy * oy + oz * oz);
            float inside = std::min(std::max({qx, qy, qz}), 0.0f);
            return outside + inside;
        }
        case SdfOp::RoundedBox: {
            float qx = std::fabs(px) - node.v0.x;
            float qy = std::fabs(py) - node.v0.y;
            float qz = std::fabs(pz) - node.v0.z;
            float ox = std::max(qx, 0.0f);
            float oy = std::max(qy, 0.0f);
            float oz = std::max(qz, 0.0f);
            float outside = std::sqrt(ox * ox + oy * oy + oz * oz);
            float inside = std::min(std::max({qx, qy, qz}), 0.0f);
            return outside + inside - node.f0;
        }
        case SdfOp::Cylinder: {
            float dx = std::sqrt(px * px + pz * pz) - node.f0;
            float dy = std::fabs(py) - node.f1;
            float ox = std::max(dx, 0.0f);
            float oy = std::max(dy, 0.0f);
            return std::sqrt(ox * ox + oy * oy) + std::min(std::max(dx, dy), 0.0f);
        }
        case SdfOp::Capsule: {
            float pax = px - node.v0.x, pay = py - node.v0.y, paz = pz - node.v0.z;
            float bax = node.v1.x - node.v0.x, bay = node.v1.y - node.v0.y, baz = node.v1.z - node.v0.z;
            float ba_len2 = bax * bax + bay * bay + baz * baz;
            float h = 0.0f;
            if (ba_len2 > 1e-12f) {
                h = std::clamp((pax * bax + pay * bay + paz * baz) / ba_len2, 0.0f, 1.0f);
            }
            float rx = pax - bax * h, ry = pay - bay * h, rz = paz - baz * h;
            return std::sqrt(rx * rx + ry * ry + rz * rz) - node.f0;
        }
        case SdfOp::Torus: {
            float qx = std::sqrt(px * px + pz * pz) - node.f0;
            float qy = py;
            return std::sqrt(qx * qx + qy * qy) - node.f1;
        }
        case SdfOp::Plane: {
            return (px * node.v0.x + py * node.v0.y + pz * node.v0.z) + node.f0;
        }
        case SdfOp::Union: {
            float d1 = evalScalarNode(graph, node.left, px, py, pz);
            float d2 = evalScalarNode(graph, node.right, px, py, pz);
            return std::min(d1, d2);
        }
        case SdfOp::Intersection: {
            float d1 = evalScalarNode(graph, node.left, px, py, pz);
            float d2 = evalScalarNode(graph, node.right, px, py, pz);
            return std::max(d1, d2);
        }
        case SdfOp::Subtraction: {
            float d1 = evalScalarNode(graph, node.left, px, py, pz);
            float d2 = evalScalarNode(graph, node.right, px, py, pz);
            return std::max(d1, -d2);
        }
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
        case SdfOp::Translate: {
            return evalScalarNode(graph, node.left, px - node.v0.x, py - node.v0.y, pz - node.v0.z);
        }
        case SdfOp::RotateY: {
            float c = std::cos(node.f0);
            float s = std::sin(node.f0);
            float rx = px * c - pz * s;
            float ry = py;
            float rz = px * s + pz * c;
            return evalScalarNode(graph, node.left, rx, ry, rz);
        }
        case SdfOp::ScaleUniform: {
            float s = node.f0;
            float dist = evalScalarNode(graph, node.left, px / s, py / s, pz / s);
            return dist * s;
        }
        case SdfOp::Noise3D: {
            float nx = px * node.v0.x;
            float ny = py * node.v0.y;
            float nz = pz * node.v0.z;
            float flrx = std::floor(nx), flry = std::floor(ny), flrz = std::floor(nz);
            int32_t i0 = static_cast<int32_t>(flrx);
            int32_t j0 = static_cast<int32_t>(flry);
            int32_t k0 = static_cast<int32_t>(flrz);
            float fx = nx - flrx, fy = ny - flry, fz = nz - flrz;
            float u = fx * fx * (3.0f - 2.0f * fx);
            float v = fy * fy * (3.0f - 2.0f * fy);
            float w = fz * fz * (3.0f - 2.0f * fz);
            auto hashc = [](int32_t ix, int32_t iy, int32_t iz) -> float {
                uint32_t h = static_cast<uint32_t>(ix) ^ (static_cast<uint32_t>(iy) * 374761393u);
                h ^= static_cast<uint32_t>(iz) * 668265263u;
                h = (h ^ (h >> 13)) * 1274126177u;
                h ^= h >> 16;
                uint32_t m = h & 0x7FFFu;
                return static_cast<float>(m) * (2.0f / 32767.0f) - 1.0f;
            };
            float h000 = hashc(i0, j0, k0), h100 = hashc(i0 + 1, j0, k0);
            float h010 = hashc(i0, j0 + 1, k0), h110 = hashc(i0 + 1, j0 + 1, k0);
            float h001 = hashc(i0, j0, k0 + 1), h101 = hashc(i0 + 1, j0, k0 + 1);
            float h011 = hashc(i0, j0 + 1, k0 + 1), h111 = hashc(i0 + 1, j0 + 1, k0 + 1);
            float x00 = h000 + (h100 - h000) * u, x10 = h010 + (h110 - h010) * u;
            float x01 = h001 + (h101 - h001) * u, x11 = h011 + (h111 - h011) * u;
            float y0 = x00 + (x10 - x00) * v, y1 = x01 + (x11 - x01) * v;
            float n = y0 + (y1 - y0) * w;
            return n * node.f0;
        }
        case SdfOp::Displace: {
            float d = evalScalarNode(graph, node.left, px, py, pz);
            float n = evalScalarNode(graph, node.right, px, py, pz);
            return d + n;
        }
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
            opts.enable_optimizations = true;
            opts.enable_avx2 = true;
            opts.enable_fma = true;

            // 1. Grid Evaluation Function
            std::string gridModName = "sdf_grid_" + std::to_string(hash);
            Module modGrid(gridModName);
            Function* fnGrid = modGrid.create_function("eval_grid", Type::void_type(), {
                Type::ptr(), // float* out_field
                Type::ptr(), // const float* bounds_min
                Type::ptr(), // const float* cell_size
                Type::i32(), // int32_t dim_x
                Type::i32(), // int32_t dim_y
                Type::i32()  // int32_t dim_z
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
                return lowerToMir(kbGrid, graph, root, px, py, pz);
            });
            kbGrid.builder().build_ret_void();

            KernelJit jitGrid(opts);
            kernel.gridKernel = jitGrid.compile(*fnGrid);
            if (kernel.gridKernel.is_valid()) {
                kernel.gridFn = kernel.gridKernel.as<SdfGridEvalFn>();
            }

            // 2. Point Evaluation Function
            std::string pointModName = "sdf_point_" + std::to_string(hash);
            Module modPoint(pointModName);
            Function* fnPoint = modPoint.create_function("eval_point", Type::f32(), {
                Type::f32(), Type::f32(), Type::f32()
            });

            SdfKernelBuilder kbPoint(modPoint, fnPoint);
            BasicBlock* entryPoint = kbPoint.builder().append_block("entry");
            Value* px = kbPoint.builder().add_block_param(entryPoint, Type::f32());
            Value* py = kbPoint.builder().add_block_param(entryPoint, Type::f32());
            Value* pz = kbPoint.builder().add_block_param(entryPoint, Type::f32());

            kbPoint.position_at_end(entryPoint);
            Value* d = lowerToMir(kbPoint, graph, root, px, py, pz);
            kbPoint.builder().build_ret(d);

            KernelJit jitPoint(opts);
            kernel.pointKernel = jitPoint.compile(*fnPoint);
            if (kernel.pointKernel.is_valid()) {
                kernel.pointFn = kernel.pointKernel.as<SdfPointEvalFn>();
            }

            kernel.isJit = (kernel.gridFn != nullptr && kernel.pointFn != nullptr);
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

SDFVolume JitSdfCompiler::evaluateVolume(const SdfGraph& graph, int dimX, int dimY, int dimZ,
                                        const bromath::AABB3& bounds, int rootNode) {
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

    if (kernel.gridFn) {
        float boundsMin[3] = {bounds.min.x, bounds.min.y, bounds.min.z};
        kernel.gridFn(vol.field.data(), boundsMin, vol.cellSize, dimX, dimY, dimZ);
    } else {
        // Fallback scalar evaluation
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
