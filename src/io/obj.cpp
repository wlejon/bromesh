#include "bromesh/io/obj.h"

namespace bromesh {

MeshData loadOBJ(const std::string& path) {
    (void)path;
    return {};
}

bool saveOBJ(const MeshData& mesh, const std::string& path) {
    (void)mesh; (void)path;
    return false;
}

} // namespace bromesh
