#include "bromesh/analysis/bbox.h"
#include <limits>

namespace bromesh {

BBox computeBBox(const MeshData& mesh) {
    (void)mesh;
    return {};
}

bool isManifold(const MeshData& mesh) {
    (void)mesh;
    return false;
}

float computeVolume(const MeshData& mesh) {
    (void)mesh;
    return 0.0f;
}

} // namespace bromesh
