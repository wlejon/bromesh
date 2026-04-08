#include "bromesh/primitives/primitives.h"

namespace bromesh {

MeshData box(float halfW, float halfH, float halfD) {
    (void)halfW; (void)halfH; (void)halfD;
    return {};
}

MeshData sphere(float radius, int segments, int rings) {
    (void)radius; (void)segments; (void)rings;
    return {};
}

MeshData cylinder(float radius, float halfHeight, int segments) {
    (void)radius; (void)halfHeight; (void)segments;
    return {};
}

MeshData capsule(float radius, float halfHeight, int segments, int rings) {
    (void)radius; (void)halfHeight; (void)segments; (void)rings;
    return {};
}

MeshData plane(float halfW, float halfD, int subdivX, int subdivZ) {
    (void)halfW; (void)halfD; (void)subdivX; (void)subdivZ;
    return {};
}

MeshData torus(float majorRadius, float minorRadius, int majorSegments, int minorSegments) {
    (void)majorRadius; (void)minorRadius; (void)majorSegments; (void)minorSegments;
    return {};
}

MeshData heightmapGrid(const float* heights, int gridW, int gridH, float cellSize) {
    (void)heights; (void)gridW; (void)gridH; (void)cellSize;
    return {};
}

} // namespace bromesh
