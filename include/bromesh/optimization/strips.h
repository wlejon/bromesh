#pragma once

#include <vector>
#include <cstdint>

namespace bromesh {

/// Convert a triangle list index buffer to a triangle strip.
/// restartIndex: use 0xFFFFFFFF for restart indices, or 0 for degenerate triangles.
/// Returns the strip index buffer (may be larger than input due to restarts/degenerates).
std::vector<uint32_t> stripify(const std::vector<uint32_t>& indices,
                               size_t vertexCount,
                               uint32_t restartIndex = 0xFFFFFFFF);

/// Convert a triangle strip index buffer back to a triangle list.
/// restartIndex: must match the value used during stripification.
std::vector<uint32_t> unstripify(const std::vector<uint32_t>& strip,
                                 uint32_t restartIndex = 0xFFFFFFFF);

} // namespace bromesh
