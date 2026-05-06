#pragma once

#include "bromesh/mesh_data.h"

#include <string>
#include <vector>

namespace bromesh {

/// Leaf shape — selects the UV cell in a 4x4 atlas grid.
/// Cell indices (column, row) with origin top-left:
///   Oval    -> (0,0)
///   Pointed -> (1,0)
///   Lobed   -> (2,0)
///   Needle  -> (3,0)
///   Frond   -> (0,1)
///   Petal   -> (1,1)
enum class LeafShape : int {
    Oval = 0,
    Pointed = 1,
    Lobed = 2,
    Needle = 3,
    Frond = 4,
    Petal = 5,
};

struct LeafCardOptions {
    /// Card width along local X (full width).
    float width = 0.4f;
    /// Card length along local Z (full length).
    float length = 1.0f;
    /// Bend amount along the length axis (radians at the tip — the card
    /// curls forward as a circular arc with this much total deflection).
    float bend = 0.0f;
    /// Twist around the length (Z) axis, radians from base to tip.
    float curl = 0.0f;
    /// If true, the card pivots at the base (z = 0). If false, pivot at center.
    bool stemOffset = true;
    /// Number of segments along width (>= 1).
    int widthSegments = 4;
    /// Number of segments along length (>= 1). Bend looks smooth at >= 8.
    int lengthSegments = 8;
    /// If true, UVs span the full [0,1] range instead of selecting a sub-cell
    /// of the 4x4 atlas. Use when the consumer binds a per-leaf texture
    /// (procedural petal/leaf with built-in alpha mask) rather than a packed
    /// atlas. The `shape` parameter is ignored when this is true.
    bool fullUV = false;
};

/// Build a low-poly leaf or petal card with UVs sampled from a 4x4 atlas.
/// Output mesh is in local space, lying in the XZ plane (Y is the card normal).
/// Vertex colors are populated with windBend in the R channel: 0 at the base,
/// 1 at the tip (consumers wire this into the wind shader).
MeshData leafCard(LeafShape shape, const LeafCardOptions& opts = {});

struct FlowerOptions {
    int petalCount = 6;
    LeafShape petalShape = LeafShape::Petal;
    float petalLength = 0.5f;
    float petalWidth = 0.25f;
    float petalCurl = 0.0f;       // radians per petal twist
    float petalBend = 0.6f;       // tip deflection radians
    int layers = 1;               // >1 stacks rings (rose-like)
    float layerTwist = 0.4f;      // radians between successive layers
    float centerRadius = 0.08f;
    float centerHeight = 0.04f;
    float centerColor[3] = {1.0f, 0.85f, 0.2f};
};

/// Build a radial flower: a small dome center with `petalCount` leaf cards
/// arranged in `layers` rings. Each subsequent ring shrinks slightly and is
/// rotated by `layerTwist`. Returns a single merged MeshData (positions,
/// normals, UVs, colors).
MeshData flower(const FlowerOptions& opts = {});

} // namespace bromesh
