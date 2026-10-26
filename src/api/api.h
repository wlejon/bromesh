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
const bromesh::MeshData* meshDataOf(bronze::Value v);
bromesh::MeshData* meshDataOfMut(bronze::Value v);
const bromesh::SkinData* skinDataOf(bronze::Value v);
const bromesh::Skeleton* skeletonOf(bronze::Value v);
const bromesh::Animation* animationOf(bronze::Value v);

/// How a path handed to the file loaders and savers (`Mesh.loadOBJ`,
/// `mesh.saveGLTF`, ...) becomes a filesystem path. Unset, the path is used
/// as given; a host sets its `fs` resolver so a relative path means what it
/// means to the app. Process-wide: every realm shares the host's filesystem.
void setPathResolver(std::function<std::string(const std::string&)> resolver);

// ---- Access -----------------------------------------------------------------
//
// Which files a program may touch is the host's decision, not bromesh's: every
// native that takes a path asks the host's one access gate,
// bronze::embed::checkAccess (embed/embed_access.h), with the path resolved
// (setPathResolver), absolute and lexically normal (UTF-8), on the JS thread
// that made the call, before the file is opened. With no gate installed every
// path is allowed and nothing below changes.
//
// What each native asks:
///   Mesh.loadOBJ / loadPLY / loadSTL / loadFBX / loadVOX / loadSplatPLY,
///   Mesh.loadGLTF (bro.mesh.loadGLTF), Rig.specFromFile (bro.rigging.specFromFile)
///                 the path Read (a glTF's buffers and images beside it are the
///                 library's to open: the file the app named is what is asked)
///   mesh.saveOBJ / savePLY / saveSTL / saveGLTF, Mesh.saveSplatPLY
///                 the path Write
//
// A refusal is bronze::embed::newAccessError, thrown: an Error whose `message`
// is the gate's sentence (named PermissionRevoked for a grant taken back),
// `code` is "EACCES" and `path` the absolute path asked about. Nothing is read
// or written.

} // namespace bromesh::api

using bromesh::api::installMesh;
using bromesh::api::installRigging;
