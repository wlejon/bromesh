#include "bromesh/optimization/strips.h"

#ifdef BROMESH_HAS_MESHOPTIMIZER
#include <meshoptimizer.h>
#endif

namespace bromesh {

std::vector<uint32_t> stripify(const std::vector<uint32_t>& indices,
                               size_t vertexCount,
                               uint32_t restartIndex) {
#ifdef BROMESH_HAS_MESHOPTIMIZER
    if (indices.empty()) return {};

    size_t bound = meshopt_stripifyBound(indices.size());
    std::vector<uint32_t> strip(bound);
    size_t stripSize = meshopt_stripify(
        strip.data(),
        indices.data(),
        indices.size(),
        vertexCount,
        restartIndex
    );
    strip.resize(stripSize);
    return strip;
#else
    (void)indices; (void)vertexCount; (void)restartIndex;
    return {};
#endif
}

std::vector<uint32_t> unstripify(const std::vector<uint32_t>& strip,
                                 uint32_t restartIndex) {
#ifdef BROMESH_HAS_MESHOPTIMIZER
    if (strip.empty()) return {};

    size_t bound = meshopt_unstripifyBound(strip.size());
    std::vector<uint32_t> indices(bound);
    size_t indexCount = meshopt_unstripify(
        indices.data(),
        strip.data(),
        strip.size(),
        restartIndex
    );
    indices.resize(indexCount);
    return indices;
#else
    (void)strip; (void)restartIndex;
    return {};
#endif
}

} // namespace bromesh
