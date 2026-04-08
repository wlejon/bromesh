#pragma once

#include "bromesh/mesh_data.h"

namespace bromesh {

/// Isotropic remeshing: improve triangle quality by making edges uniform length.
/// targetEdgeLength: desired edge length (0 = auto from average edge length).
/// iterations: number of remeshing passes (more = smoother, slower).
/// Each pass: split long edges, collapse short edges, flip to improve valence,
/// then tangential relaxation.
MeshData remeshIsotropic(const MeshData& mesh, float targetEdgeLength = 0.0f,
                         int iterations = 5);

} // namespace bromesh
