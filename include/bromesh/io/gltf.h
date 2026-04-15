#pragma once

#include "bromesh/mesh_data.h"
#include <string>
#include <vector>

namespace bromesh {

/// Result of loading a glTF/GLB file.
///
/// Parallel arrays: for each mesh i,
///   - skins[i] holds its per-vertex bone weights/indices (empty if unskinned)
///   - meshSkeleton[i] is the index into skeletons (-1 if unskinned)
///
/// Animations target skeletons by index via animationSkeleton[k]. Most files
/// have a single skeleton and all animations bound to it.
struct GltfScene {
    std::vector<MeshData> meshes;
    std::vector<SkinData> skins;
    std::vector<int>      meshSkeleton;
    std::vector<Skeleton> skeletons;
    std::vector<Animation> animations;
    std::vector<int>       animationSkeleton;
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
