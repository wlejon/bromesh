#pragma once

#include "bromesh/mesh_data.h"
#include <string>
#include <vector>

namespace bromesh {

/// Result of loading a glTF/GLB file (may contain multiple meshes).
struct GltfScene {
    std::vector<MeshData> meshes;
    std::vector<SkinData> skins;
    std::vector<MorphTarget> morphTargets;
};

/// Load a glTF or GLB file. Returns empty scene on failure.
GltfScene loadGLTF(const std::string& path);

/// Save a mesh as glTF (.gltf + .bin) or GLB (.glb, determined by extension).
bool saveGLTF(const MeshData& mesh, const std::string& path);

} // namespace bromesh
