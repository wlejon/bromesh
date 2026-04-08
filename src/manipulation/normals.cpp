#include "bromesh/manipulation/normals.h"

namespace bromesh {

void computeNormals(MeshData& mesh) {
    (void)mesh;
}

MeshData computeFlatNormals(const MeshData& mesh) {
    (void)mesh;
    return {};
}

std::vector<float> computeTangents(const MeshData& mesh) {
    (void)mesh;
    return {};
}

} // namespace bromesh
