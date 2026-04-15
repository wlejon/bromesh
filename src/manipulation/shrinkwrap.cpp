#include "bromesh/manipulation/shrinkwrap.h"
#include "bromesh/analysis/bvh.h"
#include "bromesh/analysis/raycast.h"

namespace bromesh {

static void tryProject(MeshData& mesh, size_t vi,
                       const MeshData& target, const MeshBVH& bvh,
                       const float* dir, float maxDistance, float offset) {
    float origin[3] = {
        mesh.positions[vi * 3 + 0],
        mesh.positions[vi * 3 + 1],
        mesh.positions[vi * 3 + 2],
    };

    // Try +dir then -dir (start just behind origin to handle vertices already on surface)
    float back[3] = { origin[0] - dir[0] * 1e-3f, origin[1] - dir[1] * 1e-3f, origin[2] - dir[2] * 1e-3f };
    RayHit hit = bvh.raycast(target, back, dir, maxDistance);
    if (!hit.hit) {
        float neg[3] = { -dir[0], -dir[1], -dir[2] };
        float back2[3] = { origin[0] - neg[0] * 1e-3f, origin[1] - neg[1] * 1e-3f, origin[2] - neg[2] * 1e-3f };
        hit = bvh.raycast(target, back2, neg, maxDistance);
    }
    if (!hit.hit) return;

    mesh.positions[vi * 3 + 0] = hit.position[0] + hit.normal[0] * offset;
    mesh.positions[vi * 3 + 1] = hit.position[1] + hit.normal[1] * offset;
    mesh.positions[vi * 3 + 2] = hit.position[2] + hit.normal[2] * offset;
}

void shrinkwrap(MeshData& mesh,
                const MeshData& target,
                ShrinkwrapMode mode,
                float maxDistance,
                float offset,
                const float axis[3]) {
    if (mesh.empty() || target.empty()) return;
    MeshBVH bvh = MeshBVH::build(target);

    const bool hasNormals = mesh.hasNormals();
    float axisVec[3] = { 0, 1, 0 };
    if (axis) { axisVec[0] = axis[0]; axisVec[1] = axis[1]; axisVec[2] = axis[2]; }

    for (size_t vi = 0; vi < mesh.vertexCount(); ++vi) {
        if (mode == ShrinkwrapMode::Nearest) {
            float p[3] = {
                mesh.positions[vi * 3 + 0],
                mesh.positions[vi * 3 + 1],
                mesh.positions[vi * 3 + 2],
            };
            RayHit h = closestPoint(target, p);
            if (!h.hit) continue;
            if (maxDistance > 0 && h.distance > maxDistance) continue;
            mesh.positions[vi * 3 + 0] = h.position[0] + h.normal[0] * offset;
            mesh.positions[vi * 3 + 1] = h.position[1] + h.normal[1] * offset;
            mesh.positions[vi * 3 + 2] = h.position[2] + h.normal[2] * offset;
        } else if (mode == ShrinkwrapMode::ProjectAlongNormal && hasNormals) {
            float dir[3] = {
                mesh.normals[vi * 3 + 0],
                mesh.normals[vi * 3 + 1],
                mesh.normals[vi * 3 + 2],
            };
            tryProject(mesh, vi, target, bvh, dir, maxDistance, offset);
        } else if (mode == ShrinkwrapMode::ProjectAlongAxis) {
            tryProject(mesh, vi, target, bvh, axisVec, maxDistance, offset);
        }
    }
}

} // namespace bromesh
