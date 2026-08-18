#pragma once

#include "bromesh/mesh_data.h"
#include "bromesh/procedural/plants.h"
#include <bromath/vec.h>

namespace bromesh {

/// Botanical leaf arrangement pattern along a twig or shoot axis.
enum class Phyllotaxy : int {
    Alternate = 0,       ///< Distichous / 2-ranked spiral, alternating sides along twig (oak, elm, birch, beech)
    Opposite = 1,        ///< Decussate paired leaves opposite each other with 90° node twist (maple, ash, lilac)
    Spiral = 2,          ///< Golden angle (~137.5°) rosette along shoot axis (magnolia, apple, cherry)
    Fascicle = 3,        ///< Pine needle bundle: 2–5 needles radiating from a basal sheath
    CompoundPinnate = 4, ///< Paired lateral leaflets along a central rachis with a terminal leaflet (walnut, rowan, acacia)
};

/// Configuration options for procedural botanical leaf clusters and twig sprays.
struct LeafClusterOptions {
    /// Number of leaves / leaflets in the cluster.
    int count = 6;
    /// Length of the supporting micro-twig / rachis along local +Z.
    float twigLength = 0.25f;
    /// Radius of the supporting micro-twig.
    float twigRadius = 0.005f;
    /// Length of the individual leaf stalk (petiole) connecting leaf to twig.
    float petioleLength = 0.04f;
    /// Width of individual leaf cards.
    float leafWidth = 0.12f;
    /// Length of individual leaf cards.
    float leafLength = 0.20f;
    /// Leaf shape profile (atlas cell or silhouette).
    LeafShape shape = LeafShape::Oval;
    /// Length-wise bend deflection (radians).
    float leafBend = 0.3f;
    /// Axial twist curl (radians).
    float leafCurl = 0.1f;
    /// Bilateral transverse cupping.
    float leafCup = 0.2f;
    /// Gravitational sag deflection along petiole and leaf.
    float droop = 0.2f;
    /// Phototropic bias: adaxial surface turned toward sky/light (+Y).
    float upBias = 0.6f;
    /// Lateral fan/divergence angle (radians) from the central twig axis.
    float spread = 0.7f;
    /// If true, generates the micro-twig cylinder stem so leaves aren't floating in space.
    bool includeTwigMesh = true;
    /// If true, modulates leaf card width geometrically by shape silhouette.
    bool shapedSilhouette = true;
    /// If true, leaf card UVs span full [0, 1] instead of 4x4 atlas cell.
    bool fullUV = false;
};

/// Build a low-poly botanical leaf cluster / twig spray with petioles and leaves
/// in the specified phyllotaxy arrangement.
///
/// Output mesh is in local space:
/// - Twig root is at (0, 0, 0), extending along local +Z to (0, 0, twigLength)
/// - Local +Y is the upward / light-facing normal direction
/// - Local +X is the lateral spreading axis
///
/// Vertex colors encode wind bend in the R channel: 0.0 at the twig base attachment,
/// scaling smoothly to 1.0 at the outermost leaf tips.
MeshData leafCluster(Phyllotaxy phyllotaxy, const LeafClusterOptions& opts = {});

} // namespace bromesh
