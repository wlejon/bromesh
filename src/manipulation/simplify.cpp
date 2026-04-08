#include "bromesh/manipulation/simplify.h"

namespace bromesh {

MeshData simplify(const MeshData& mesh, float targetRatio, float targetError) {
    (void)mesh; (void)targetRatio; (void)targetError;
    return {};
}

std::vector<MeshData> generateLODChain(const MeshData& mesh,
                                       const float* ratios, int count) {
    (void)mesh; (void)ratios; (void)count;
    return {};
}

} // namespace bromesh
