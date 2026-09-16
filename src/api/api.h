#pragma once

#include "embed/embed.h"

namespace bromesh::api {

/// Mounts `Mesh`, `MeshBVH`, `ProgressiveMesh`, and `bro.mesh` in the current Bronze realm.
void installMesh();

/// Mounts `Skeleton`, `Joint`, `SkeletonRig`, `AnimationClip`, `Pose`, `SkinData`, `IK`,
/// and `bro.rigging` in the current Bronze realm.
void installRigging();

} // namespace bromesh::api

using bromesh::api::installMesh;
using bromesh::api::installRigging;
