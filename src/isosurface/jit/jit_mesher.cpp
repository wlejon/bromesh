#include "bromesh/isosurface/jit/jit_mesher.h"
#include "bromesh/isosurface/jit/sdf_compiler.h"
#include "bromesh/isosurface/marching_cubes.h"
#include "bromesh/isosurface/surface_nets.h"
#include "bromesh/manipulation/normals.h"

#include <algorithm>
#include <cmath>

namespace bromesh {

static void computeAnalyticalNormals(MeshData& mesh, const CompiledSdfKernel& kernel, const float* cellSize) {
    if ((!kernel.gradFn && !kernel.pointFn) || mesh.positions.empty()) {
        computeNormals(mesh);
        return;
    }

    mesh.normals.resize(mesh.positions.size());
    const float minCell = std::min({cellSize[0], cellSize[1], cellSize[2]});
    const float eps = std::max(1e-4f, 0.25f * minCell);
    const float inv2eps = 1.0f / (2.0f * eps);
    const int64_t numVerts = static_cast<int64_t>(mesh.positions.size() / 3);
    SdfGradEvalFn gradFn = kernel.gradFn;
    SdfPointEvalFn pointFn = kernel.pointFn;

    auto evalVertRange = [&](int64_t start, int64_t end) {
        if (gradFn) {
            for (int64_t i = start; i < end; ++i) {
                float x = mesh.positions[3 * i + 0];
                float y = mesh.positions[3 * i + 1];
                float z = mesh.positions[3 * i + 2];

                float gx = 0.0f, gy = 0.0f, gz = 0.0f;
                gradFn(x, y, z, &gx, &gy, &gz);

                float len = std::sqrt(gx * gx + gy * gy + gz * gz);
                if (len > 1e-7f) {
                    float invLen = 1.0f / len;
                    mesh.normals[3 * i + 0] = gx * invLen;
                    mesh.normals[3 * i + 1] = gy * invLen;
                    mesh.normals[3 * i + 2] = gz * invLen;
                } else {
                    mesh.normals[3 * i + 0] = 0.0f;
                    mesh.normals[3 * i + 1] = 1.0f;
                    mesh.normals[3 * i + 2] = 0.0f;
                }
            }
        } else {
            for (int64_t i = start; i < end; ++i) {
                float x = mesh.positions[3 * i + 0];
                float y = mesh.positions[3 * i + 1];
                float z = mesh.positions[3 * i + 2];

                float gx = (pointFn(x + eps, y, z) - pointFn(x - eps, y, z)) * inv2eps;
                float gy = (pointFn(x, y + eps, z) - pointFn(x, y - eps, z)) * inv2eps;
                float gz = (pointFn(x, y, z + eps) - pointFn(x, y, z - eps)) * inv2eps;

                float len = std::sqrt(gx * gx + gy * gy + gz * gz);
                if (len > 1e-7f) {
                    float invLen = 1.0f / len;
                    mesh.normals[3 * i + 0] = gx * invLen;
                    mesh.normals[3 * i + 1] = gy * invLen;
                    mesh.normals[3 * i + 2] = gz * invLen;
                } else {
                    mesh.normals[3 * i + 0] = 0.0f;
                    mesh.normals[3 * i + 1] = 1.0f;
                    mesh.normals[3 * i + 2] = 0.0f;
                }
            }
        }
    };

#if BROMESH_HAS_BRASS_JIT
    if (numVerts >= 64) {
        brass::codegen::parallel_for(0, numVerts, evalVertRange);
    } else {
        evalVertRange(0, numVerts);
    }
#else
    evalVertRange(0, numVerts);
#endif
}

MeshData marchingCubesFromSDF(
    const SdfGraph& graph,
    int dimX, int dimY, int dimZ,
    const bromath::AABB3& bounds,
    float isoLevel,
    bool closeBoundary,
    bool computeGradients,
    bool enableMacroBlockPruning
) {
    SDFVolume vol = JitSdfCompiler::instance().evaluateVolume(
        graph, dimX, dimY, dimZ, bounds, -1, true, isoLevel, enableMacroBlockPruning);
    if (vol.field.empty()) {
        return MeshData();
    }

    const uint8_t* blockMask = (!vol.activeBlockMask.empty() && vol.activeBlockCount < vol.totalBlockCount)
                                   ? vol.activeBlockMask.data()
                                   : nullptr;

    MeshData mesh = marchingCubes(vol.field.data(), dimX, dimY, dimZ, isoLevel, 1.0f, closeBoundary, blockMask);
    if (mesh.empty()) {
        return mesh;
    }

    // Transform vertices from grid-index coordinate frame to world-space bounds
    const size_t numVerts = mesh.positions.size() / 3;
    for (size_t i = 0; i < numVerts; ++i) {
        mesh.positions[3 * i + 0] = bounds.min.x + mesh.positions[3 * i + 0] * vol.cellSize[0];
        mesh.positions[3 * i + 1] = bounds.min.y + mesh.positions[3 * i + 1] * vol.cellSize[1];
        mesh.positions[3 * i + 2] = bounds.min.z + mesh.positions[3 * i + 2] * vol.cellSize[2];
    }

    if (computeGradients) {
        CompiledSdfKernel kernel = JitSdfCompiler::instance().compile(graph);
        computeAnalyticalNormals(mesh, kernel, vol.cellSize);
    } else {
        computeNormals(mesh);
    }

    return mesh;
}

MeshData surfaceNetsFromSDF(
    const SdfGraph& graph,
    int dimX, int dimY, int dimZ,
    const bromath::AABB3& bounds,
    float isoLevel,
    bool enableMacroBlockPruning
) {
    SDFVolume vol = JitSdfCompiler::instance().evaluateVolume(
        graph, dimX, dimY, dimZ, bounds, -1, true, isoLevel, enableMacroBlockPruning);
    if (vol.field.empty()) {
        return MeshData();
    }

    const uint8_t* blockMask = (!vol.activeBlockMask.empty() && vol.activeBlockCount < vol.totalBlockCount)
                                   ? vol.activeBlockMask.data()
                                   : nullptr;

    MeshData mesh = surfaceNets(vol.field.data(), dimX, dimY, dimZ, isoLevel, 1.0f, blockMask);
    if (mesh.empty()) {
        return mesh;
    }

    // Transform vertices from grid-index coordinate frame to world-space bounds
    const size_t numVerts = mesh.positions.size() / 3;
    for (size_t i = 0; i < numVerts; ++i) {
        mesh.positions[3 * i + 0] = bounds.min.x + mesh.positions[3 * i + 0] * vol.cellSize[0];
        mesh.positions[3 * i + 1] = bounds.min.y + mesh.positions[3 * i + 1] * vol.cellSize[1];
        mesh.positions[3 * i + 2] = bounds.min.z + mesh.positions[3 * i + 2] * vol.cellSize[2];
    }

    CompiledSdfKernel kernel = JitSdfCompiler::instance().compile(graph);
    computeAnalyticalNormals(mesh, kernel, vol.cellSize);

    return mesh;
}

} // namespace bromesh
