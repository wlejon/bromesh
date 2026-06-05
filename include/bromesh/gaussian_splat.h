#pragma once

#include <bromath/aabb.h>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace bromesh {

/// A 3D Gaussian Splat cloud — the output of 3DGS reconstruction / generation
/// (e.g. TripoSplat) and the input to the scene GaussianSplatNode rasterizer.
///
/// Separate attribute streams (SoA), mirroring MeshData, so each maps to a JS
/// TypedArray with zero repacking.
///
/// Values are stored RENDER-READY (activated), not in the raw .ply encoding:
///   - scales    : linear world-space std-devs along the local axes (the .ply's
///                 log-scale already exponentiated).
///   - opacities : [0,1] alpha (the .ply's logit already passed through sigmoid).
///   - rotations : unit quaternions xyzw (normalized).
///   - sh        : spherical-harmonic color coefficients, left as coefficients
///                 (that *is* the render-ready form — the rasterizer evaluates
///                 SH -> RGB per view direction).
/// loadSplatPLY / saveSplatPLY convert to/from the raw log-scale / logit-opacity
/// storage that standard 3DGS .ply files (INRIA / PlayCanvas) use.
///
/// SH memory layout (`sh`): interleaved by coefficient, RGB together —
///   [r0 g0 b0  r1 g1 b1  ...  r_{C-1} g_{C-1} b_{C-1}]   per splat,
/// where coefficient 0 is the DC term (base color) and C = (shDegree+1)^2.
/// stride = shStride() = 3 * C. This differs from the .ply on-disk order, which
/// is channel-major (all of R's coeffs, then G's, then B's); the loader/saver
/// transpose between the two.
struct GaussianSplatCloud {
    std::vector<float> positions; // xyz,  stride 3
    std::vector<float> scales;    // xyz,  stride 3 (linear std-dev)
    std::vector<float> rotations; // xyzw, stride 4 (unit quaternion)
    std::vector<float> opacities; // a,    stride 1 ([0,1])
    std::vector<float> sh;        // stride shStride(); see layout note above
    int shDegree = 0;             // 0..3

    size_t count() const { return positions.size() / 3; }
    int shStride() const { return 3 * (shDegree + 1) * (shDegree + 1); }
    bool empty() const { return positions.empty(); }

    void clear() {
        positions.clear();
        scales.clear();
        rotations.clear();
        opacities.clear();
        sh.clear();
        shDegree = 0;
    }

    void reserve(size_t n) {
        positions.reserve(n * 3);
        scales.reserve(n * 3);
        rotations.reserve(n * 4);
        opacities.reserve(n);
        sh.reserve(n * static_cast<size_t>(shStride()));
    }

    /// Axis-aligned bounds of all splat centers (centers only — does not inflate
    /// by per-splat scale). Empty/inverted box when the cloud is empty.
    bromath::AABB3 bounds() const;
};

} // namespace bromesh
