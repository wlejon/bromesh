#pragma once
#include "bromesh/mesh_data.h"
#include <bromath/vec.h>
#include <vector>

namespace bromesh {

struct LatheOptions {
    int segments = 32;               // radial angular subdivisions (>= 3)
    float startAngle = 0.0f;         // radians
    float endAngle = 6.28318530718f; // 2 * pi for full revolution
    bool capStart = true;            // cap at startAngle when partial revolution (< 2*pi)
    bool capEnd = true;              // cap at endAngle when partial revolution (< 2*pi)
    int axis = 1;                    // 0=X, 1=Y, 2=Z revolution axis
    bool smoothNormals = true;
};

/// Revolve a 2D profile curve around an axis to produce a surface of revolution
/// (vases, bottles, columns, goblets, wheels, barrels, domes).
/// Profile is a polyline in the revolution plane (x = distance from axis, y = height along axis).
MeshData lathe(const std::vector<bromath::Vec2>& profile, const LatheOptions& opts = {});

} // namespace bromesh
