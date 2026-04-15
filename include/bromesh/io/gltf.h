#pragma once

#include "bromesh/mesh_data.h"
#include <cstdint>
#include <string>
#include <vector>

namespace bromesh {

/// Decoded image payload from a glTF. RGBA8, top-left origin, tightly packed
/// (stride = width * 4).
struct Image {
    std::string name;               ///< original image name or "" if unnamed
    std::string mimeType;            ///< e.g. "image/png", "image/jpeg"
    int width  = 0;
    int height = 0;
    std::vector<uint8_t> data;       ///< RGBA8, size = width * height * 4
};

/// Flat PBR material subset — enough for the baseColor texture + tint that
/// MeshyAI and similar exports produce. Indexes refer to GltfScene::images.
struct Material {
    std::string name;
    float baseColorFactor[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    int   baseColorTexture   = -1;   ///< index into GltfScene::images, -1 if none
};

/// Result of loading a glTF/GLB file.
///
/// Parallel arrays: for each mesh i,
///   - skins[i] holds its per-vertex bone weights/indices (empty if unskinned)
///   - meshSkeleton[i] is the index into skeletons (-1 if unskinned)
///   - meshMaterial[i] is the index into materials (-1 if untextured)
///
/// Animations target skeletons by index via animationSkeleton[k]. Most files
/// have a single skeleton and all animations bound to it.
struct GltfScene {
    std::vector<MeshData>  meshes;
    std::vector<SkinData>  skins;
    std::vector<int>       meshSkeleton;
    std::vector<int>       meshMaterial;
    std::vector<Skeleton>  skeletons;
    std::vector<Animation> animations;
    std::vector<int>       animationSkeleton;
    std::vector<Material>  materials;
    std::vector<Image>     images;
};

/// Load a glTF or GLB file. Returns empty scene on failure.
/// Populates meshes, skins, skeletons, and animations. Each mesh is paired
/// with its skin (if any) at the same index; unskinned meshes get an empty
/// SkinData and meshSkeleton[i] == -1.
GltfScene loadGLTF(const std::string& path);

/// Save a single unskinned mesh. Convenience wrapper.
bool saveGLTF(const MeshData& mesh, const std::string& path);

/// Save a mesh with optional skin + skeleton + animations.
/// If skin or skeleton is null, the output is unskinned. Animations are only
/// emitted when a skeleton is provided.
bool saveGLTF(const MeshData& mesh,
              const SkinData* skin,
              const Skeleton* skeleton,
              const std::vector<Animation>& animations,
              const std::string& path);

} // namespace bromesh
