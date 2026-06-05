#pragma once

#include "bromesh/gaussian_splat.h"
#include <string>

namespace bromesh {

/// Load a 3D Gaussian Splat cloud from a .ply file (INRIA / PlayCanvas layout:
/// x,y,z, optional nx,ny,nz, f_dc_0..2, f_rest_0..N, opacity, scale_0..2,
/// rot_0..3; binary little-endian or ASCII). The SH degree is inferred from the
/// f_rest_* property count. Values are activated on load (exp scale, sigmoid
/// opacity, normalized quat) — see GaussianSplatCloud. Returns an empty cloud
/// on failure (missing file, no splat properties, unsupported format).
GaussianSplatCloud loadSplatPLY(const std::string& path);

/// Save a Gaussian Splat cloud as a standard 3DGS binary-little-endian .ply,
/// inverting the load activations (log scale, logit opacity) so the result
/// round-trips through standard 3DGS tooling and back through loadSplatPLY.
bool saveSplatPLY(const GaussianSplatCloud& cloud, const std::string& path);

} // namespace bromesh
