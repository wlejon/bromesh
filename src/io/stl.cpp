#include "bromesh/io/stl.h"

namespace bromesh {

bool saveSTL(const MeshData& mesh, const std::string& path) {
    (void)mesh; (void)path;
    return false;
}

MeshData loadSTL(const std::string& path) {
    (void)path;
    return {};
}

} // namespace bromesh
