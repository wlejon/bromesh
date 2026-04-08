#pragma once

#include "bromesh/mesh_data.h"

#include <cstdint>
#include <vector>

namespace bromesh {

/// Result of a single ray-mesh intersection.
struct RayHit {
    bool hit = false;
    float distance = 0.0f;         // distance along ray to hit point
    float position[3] = {0,0,0};   // world-space hit point
    float normal[3] = {0,0,0};     // interpolated normal at hit (if mesh has normals), else face normal
    float baryU = 0, baryV = 0, baryW = 0; // barycentric coordinates
    uint32_t triangleIndex = 0;    // which triangle was hit
};

/// Cast a ray against a mesh and return the closest hit.
/// origin: ray origin (3 floats). direction: ray direction (3 floats, need not be normalized).
/// maxDistance: maximum ray distance (0 = unlimited).
RayHit raycast(const MeshData& mesh,
               const float* origin, const float* direction,
               float maxDistance = 0.0f);

/// Cast a ray and return ALL hits (sorted by distance, nearest first).
/// Useful for transparency, thickness queries, etc.
std::vector<RayHit> raycastAll(const MeshData& mesh,
                                const float* origin, const float* direction,
                                float maxDistance = 0.0f);

/// Test if a ray hits the mesh at all (early-out, faster than raycast).
bool raycastTest(const MeshData& mesh,
                 const float* origin, const float* direction,
                 float maxDistance = 0.0f);

/// Find the closest point on the mesh surface to a query point.
/// Returns a RayHit-like result with position, normal, triangle, and barycentrics.
/// distance is the Euclidean distance from the query point to the closest surface point.
RayHit closestPoint(const MeshData& mesh, const float* point);

} // namespace bromesh
