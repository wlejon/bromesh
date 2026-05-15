#pragma once

#include "bromesh/procedural/lsystem.h"
#include "bromesh/procedural/space_colonization.h"
#include <bromath/vec.h>

#include <vector>

namespace bromesh {

/// Initial pose and per-symbol defaults for the L-system turtle.
///
/// Recognised symbols:
///   F(len?, r?)  forward; emits a segment from current pos to pos+heading*len.
///                len defaults to `stepLength`; if a second param is
///                supplied it overrides the segment radius for THIS step
///                only (the active radius is unchanged).
///   G(len?), f(len?)
///                forward without emitting a segment (move-only).
///   +(deg?)      yaw +angle around `up`
///   -(deg?)      yaw -angle
///   &(deg?)      pitch +angle around the left axis (heading × up)
///   ^(deg?)      pitch -angle
///   \\(deg?)     roll +angle around `heading`
///   /(deg?)      roll -angle
///   |            turn 180° around `up`
///   !(r)         set the active radius for subsequent F segments
///   [            push pose
///   ]            pop pose
///
/// Numeric parameters on rotation symbols are interpreted as degrees, the
/// natural unit for hand-written L-system grammars. The `angle` option in
/// `TurtleOptions` is in radians (matches the rest of bromesh) and is used
/// when a rotation symbol has no parameter.
struct TurtleOptions {
    /// Default forward step length when an `F`/`G`/`f` module has no `len`.
    float stepLength = 0.1f;
    /// Default rotation when a rotation symbol has no parameter (radians).
    float angle      = 0.4363f;     // 25 degrees
    /// Default tip radius assigned to every emitted segment unless
    /// overridden by `!` or by an `F` second parameter.
    float radius     = 0.01f;
    /// Initial pose.
    bromath::Vec3  position   = {0.0f, 0.0f, 0.0f};
    bromath::Vec3  heading    = {0.0f, 1.0f, 0.0f};   // forward
    bromath::Vec3  up         = {0.0f, 0.0f, 1.0f};   // roll reference
};

/// Interpret a flat L-system module sequence as a 3D turtle. Returns the
/// branch segments traced out, in emission order, with parent indices that
/// form a forest (one synthetic root per call when the first emitted F has
/// no parent — matches `spaceColonize`'s convention so the result drops
/// straight into `meshBranches` / `scatterLeaves`).
///
/// Unknown symbols are silently ignored, so user grammars can carry their
/// own annotations without breaking the turtle.
std::vector<BranchSegment> lsystemToBranches(
    const std::vector<Module>& modules,
    const TurtleOptions& opts = {});

} // namespace bromesh
