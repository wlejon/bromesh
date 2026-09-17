#pragma once

#include "embed/embed.h"

#include <bromesh/mesh_data.h>

#include <functional>
#include <string>

namespace bromesh::api {

/// Mounts `Mesh`, `MeshBVH`, `ProgressiveMesh`, and `bro.mesh` in the current Bronze realm.
void installMesh();

/// Mounts `Skeleton`, `Joint`, `SkeletonRig`, `AnimationClip`, `Pose`, `SkinData`, `IK`,
/// and `bro.rigging` in the current Bronze realm.
void installRigging();

/// The host's transfer seam — a Worker's postMessage moving a Mesh between
/// realms by pointer rather than by copy. Both realms install their own Mesh
/// class (one per thread, host_class.h), so the host carries the bare
/// MeshData across and asks the receiving thread for a handle of ITS class.
///
/// isMeshValue: whether `v` is a Mesh handle of the calling thread's class.
/// takeMeshData: moves the MeshData out of that handle into `out` and leaves
/// the handle holding an empty mesh (neutered, the way a transferred
/// ArrayBuffer is detached); false when `v` is not a Mesh.
/// makeMeshValue: a Mesh handle of the calling thread's class over `mesh`.
bool isMeshValue(bronze::Value v);
bool takeMeshData(bronze::Value v, bromesh::MeshData& out);
bronze::Value makeMeshValue(bromesh::MeshData mesh);

/// How a path handed to the file loaders and savers (`Mesh.loadOBJ`,
/// `mesh.saveGLTF`, ...) becomes a filesystem path. Unset, the path is used
/// as given; a host sets its `fs` resolver so a relative path means what it
/// means to the app. Process-wide: every realm shares the host's filesystem.
void setPathResolver(std::function<std::string(const std::string&)> resolver);

} // namespace bromesh::api

using bromesh::api::installMesh;
using bromesh::api::installRigging;
