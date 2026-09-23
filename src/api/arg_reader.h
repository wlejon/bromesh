#pragma once

#include "embed/embed.h"
#include "object_builder.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <span>
#include <string_view>
#include <string>
#include <vector>

namespace bromesh::api {

namespace ev = bronze::embed;
using Value = bronze::Value;

inline double numAt(std::span<const Value> args, size_t i) {
    if (i >= args.size()) return 0.0;
    Value v = args[i];
    if (ev::isObject(v)) return 0.0;
    double d = ev::toDouble(v);
    return std::isnan(d) ? 0.0 : d;
}

// The plain integer readers saturate instead of converting out of range
// (undefined behaviour for a double): 1e30 reads as INT32_MAX, not garbage.
// Counts do not use them; see intValue below.
inline int64_t i64At(std::span<const Value> args, size_t i) {
    const double d = numAt(args, i);
    if (d >= 9223372036854775807.0) return INT64_MAX;
    if (d <= -9223372036854775808.0) return INT64_MIN;
    return static_cast<int64_t>(d);
}

inline int32_t i32At(std::span<const Value> args, size_t i) {
    const double d = numAt(args, i);
    if (d >= 2147483647.0) return INT32_MAX;
    if (d <= -2147483648.0) return INT32_MIN;
    return static_cast<int32_t>(d);
}

// A double narrowed to int, saturating, NaN reading as 0. For ints where any
// value is meaningful (a tag, a parent index, a seed); counts use intValue.
inline int satInt(double d) {
    if (std::isnan(d)) return 0;
    if (d >= 2147483647.0) return INT32_MAX;
    if (d <= -2147483648.0) return INT32_MIN;
    return static_cast<int>(d);
}

inline uint32_t u32At(std::span<const Value> args, size_t i) {
    return static_cast<uint32_t>(i64At(args, i));
}

inline uint64_t u64At(std::span<const Value> args, size_t i) {
    if (i >= args.size()) return 0;
    Value v = args[i];
    if (ev::isObject(v)) return 0;
    return ev::toUint64(v);
}

// ── Validated integers ─────────────────────────────────────────────────
//
// Counts, sizes, segments, resolutions and iteration counts arrive as JS
// doubles. i32At/u32At above are for values where any int will do; a
// negative, NaN or huge double sent through them to a count wraps (a
// static_cast<size_t>(-1.0) is undefined behaviour and in practice ~1.8e19),
// which then sizes an allocation. Counts go through intValue instead: a
// non-number is a TypeError; NaN, a fraction or a value outside [lo, hi] a
// RangeError. It returns false once it has raised, and the caller returns
// ev::undefined() (the throw helpers' own result) so the exception reaches
// JS. `what` names the argument, e.g. "Mesh.sphere: segments".

inline constexpr double kMaxInt32 = 2147483647.0;
inline constexpr double kMaxUint32 = 4294967295.0;
// Ceiling for a count that sizes an output directly (points, instances,
// rays, iterations): 2^24. A bad_alloc cannot cross into compiled JS, so
// "as large as memory allows" is a crash, not an error.
inline constexpr double kMaxElements = 16777216.0;
// Ceiling for one axis of a grid (segments, rings, resolution, voxel
// dimensions): 2^12, whose square is kMaxElements.
inline constexpr double kMaxAxis = 4096.0;
// Ceiling for an iteration count (smoothing, remeshing, solver passes): 2^16.
inline constexpr double kMaxIterations = 65536.0;
// Ceiling for a recursive subdivision level, each of which multiplies the
// face count by about 4: 4^8 = 65536x the input.
inline constexpr double kMaxSubdivisions = 8.0;
// Ceiling for one axis of a 3D volume a call allocates as axis^3: 2^9.
inline constexpr double kMaxVolumeAxis = 512.0;
// Ceiling for the cell count of a 3D volume given per axis: 2^27.
inline constexpr double kMaxVolumeCells = 134217728.0;

inline std::string numberText(double d) {
    if (std::isnan(d)) return "NaN";
    if (std::isinf(d)) return d > 0 ? "Infinity" : "-Infinity";
    if (d == std::floor(d) && std::fabs(d) < 1e15) return std::to_string(static_cast<long long>(d));
    std::string s = std::to_string(d);
    while (!s.empty() && s.back() == '0') s.pop_back();
    return s;
}

inline bool intValue(Value v, std::string_view what, double lo, double hi, int64_t& out) {
    if (!ev::isNumber(v)) {
        ev::throwTypeError(std::string(what) + " must be a number");
        return false;
    }
    const double d = ev::toDouble(v);
    if (std::isnan(d) || d != std::floor(d) || d < lo || d > hi) {
        ev::throwRangeError(std::string(what) + " must be an integer in [" + numberText(lo) + ", " +
                            numberText(hi) + "], got " + numberText(d));
        return false;
    }
    out = static_cast<int64_t>(d);
    return true;
}

// An optional positional count: absent or undefined leaves `out` (the
// default) alone.
template <typename T>
inline bool countArg(std::span<const Value> args, size_t i, std::string_view what,
                     double lo, double hi, T& out) {
    if (i >= args.size() || ev::isUndefined(args[i])) return true;
    int64_t v = 0;
    if (!intValue(args[i], what, lo, hi, v)) return false;
    out = static_cast<T>(v);
    return true;
}

// A count read from an option object: an undefined (or non-object
// receiver) leaves `out` alone.
template <typename T>
inline bool countField(Value obj, std::string_view key, std::string_view what,
                       double lo, double hi, T& out) {
    if (!ev::isObject(obj)) return true;
    Value v = ev::getProperty(obj, key);
    if (ev::isUndefined(v)) return true;
    int64_t n = 0;
    if (!intValue(v, what, lo, hi, n)) return false;
    out = static_cast<T>(n);
    return true;
}

// ── Seeds ──────────────────────────────────────────────────────────────
//
// A seed is an integer in [0, hi]: hi is 2^53 - 1 (the largest integer a JS
// number holds exactly) for a 64-bit seed and the type's maximum for a
// narrower one. Undefined keeps the default; a non-number is a TypeError; a
// negative, fractional, NaN or infinite one a RangeError. (Negative seeds
// used to wrap as two's complement, and a 32-bit one modulo 2^32, so two
// different arguments silently seeded alike.)
inline constexpr double kMaxSeed64 = 9007199254740991.0;

template <typename T>
inline bool seedArg(std::span<const Value> args, size_t i, std::string_view what, double hi, T& out) {
    return countArg(args, i, what, 0.0, hi, out);
}

template <typename T>
inline bool seedField(Value obj, std::string_view key, std::string_view what, double hi, T& out) {
    return countField(obj, key, what, 0.0, hi, out);
}

// A 3D grid's cell count x * y * z against kMaxVolumeCells, each axis
// already validated.
inline bool volumeCellsOk(std::string_view what, double x, double y, double z) {
    const double cells = x * y * z;
    if (cells > kMaxVolumeCells) {
        ev::throwRangeError(std::string(what) + ": the grid is " + numberText(cells) +
                            " cells, over the " + numberText(kMaxVolumeCells) + " limit");
        return false;
    }
    return true;
}

// The `length` of an array or array-like the binding copies from: ToLength
// semantics (NaN or negative is 0), capped at 2^32 - 1, so no conversion is
// undefined and no -1 becomes a 1.8e19-element loop.
inline size_t lengthValue(Value lenV) {
    if (!ev::isNumber(lenV)) return 0;
    const double d = ev::toDouble(lenV);
    if (!(d > 0.0)) return 0;
    if (d >= kMaxUint32) return static_cast<size_t>(kMaxUint32);
    return static_cast<size_t>(d);
}

// The length of a list of records (bones, channels, sockets) the binding
// builds one entry per element for: over kMaxElements is a RangeError rather
// than a loop that allocates until it fails.
inline bool listLength(Value lenV, std::string_view what, size_t& n) {
    n = lengthValue(lenV);
    if (static_cast<double>(n) > kMaxElements) {
        ev::throwRangeError(std::string(what) + " has " + std::to_string(n) + " entries, over the " +
                            numberText(kMaxElements) + " limit");
        return false;
    }
    return true;
}

// Most elements a copy loop reserves up front; a longer (or lying) array-like
// grows the vector as it is read.
inline constexpr size_t kReserveCap = size_t{1} << 20;

// For a list reader with no error channel: whether a list of `n` elements
// may be copied. Past kMaxElements it is refused (refuseList, raised as a
// RangeError when the binding returns) and the reader answers empty, so an
// array-like claiming 2^32 - 1 elements is neither a 16 GB vector nor a
// four-billion-read loop.
inline bool copyLengthOk(size_t n) {
    if (static_cast<double>(n) <= kMaxElements) return true;
    refuseList("a list of " + std::to_string(n) + " elements is over the " +
               numberText(kMaxElements) + "-element limit");
    return false;
}

// JS ToUint32 (what storing into a Uint32Array does): NaN and infinities are
// 0, anything else truncates and wraps modulo 2^32. The narrowing then is
// defined where a static_cast of the double is not.
inline uint32_t toUint32Wrap(double d) {
    if (!std::isfinite(d)) return 0;
    double m = std::fmod(std::trunc(d), 4294967296.0);
    if (m < 0) m += 4294967296.0;
    return static_cast<uint32_t>(m);
}

inline bool boolAt(std::span<const Value> args, size_t i) {
    if (i >= args.size()) return false;
    return ev::toBool(args[i]);
}

inline std::string strAt(std::span<const Value> args, size_t i) {
    if (i >= args.size() || ev::isUndefined(args[i]) || ev::isNull(args[i]) || ev::isSymbol(args[i])) return "";
    return ev::toUtf8(args[i]);
}

inline Value argAt(std::span<const Value> args, size_t i) {
    return i < args.size() ? args[i] : ev::undefined();
}

inline bool hasArg(std::span<const Value> args, size_t i) {
    return i < args.size() && !ev::isUndefined(args[i]);
}

class ArgReader {
public:
    explicit ArgReader(std::span<const Value> args) : args_(args) {}

    double getDouble(size_t i, double def = 0.0) const {
        return hasArg(args_, i) ? numAt(args_, i) : def;
    }
    int getInt(size_t i, int def = 0) const {
        return hasArg(args_, i) ? i32At(args_, i) : def;
    }
    uint32_t getUint(size_t i, uint32_t def = 0) const {
        return hasArg(args_, i) ? u32At(args_, i) : def;
    }
    bool getBool(size_t i, bool def = false) const {
        return hasArg(args_, i) ? boolAt(args_, i) : def;
    }
    std::string getString(size_t i, const std::string& def = "") const {
        return hasArg(args_, i) ? strAt(args_, i) : def;
    }
    Value get(size_t i) const {
        return argAt(args_, i);
    }
    bool has(size_t i) const {
        return hasArg(args_, i);
    }
    size_t count() const {
        return args_.size();
    }

private:
    std::span<const Value> args_;
};

/// A GC root that reads as its CURRENT Value wherever a Value is expected.
/// For an option object read field by field: `Rooted o(a[1]);
/// objNum(o, "x", 0)` re-reads the root at every use, where `Value o = a[1]`
/// would go stale at the first allocating read.
struct Rooted {
    ev::Persistent p;
    explicit Rooted(Value v) : p(v) {}
    operator Value() const { return p.get(); }
    Value get() const { return p.get(); }
};

/// `self[name](...args)` for alias methods (computeBBox -> bounds, ...).
/// `self` is a plain copy current only at entry and the method lookup
/// allocates, so the receiver and every argument are rooted first; a throw
/// from the target propagates instead of coming back as a return value.
inline Value callMethod(Value self, std::string_view name, std::span<const Value> args) {
    ev::Persistent recv(self);
    std::vector<ev::Persistent> rooted;
    rooted.reserve(args.size());
    for (Value v : args) rooted.emplace_back(v);
    ev::Persistent fn(ev::getProperty(recv.get(), name));
    if (!ev::isFunction(fn.get())) {
        return ev::throwTypeError(std::string(name) + " is not a function");
    }
    std::vector<Value> argv;
    argv.reserve(rooted.size());
    for (const auto& p : rooted) argv.push_back(p.get());
    ev::CallResult r = ev::call(fn.get(), recv.get(), argv);
    return r.thrown ? ev::throwValue(r.value) : r.value;
}

/// Build an array from `make(i)`. The array is rooted across each make()
/// call, which usually allocates (a wrapped mesh, a nested array, a string).
inline Value hostArrayOf(size_t count, const std::function<Value(size_t)>& make) {
    ev::Persistent arr(ev::parseJson("[]").value);
    for (size_t i = 0; i < count; ++i) {
        ev::Persistent item(make(i));
        arr.set(ev::setElement(arr.get(), static_cast<uint32_t>(i), item.get()));
    }
    return arr.get();
}

/// Only for immediates (numbers, bools, undefined) or an empty span: a heap
/// Value in `items` goes stale at the first setElement, since the span is not
/// a GC root. Build heap elements inside a make() callback instead.
inline Value hostArrayOf(std::span<const Value> items) {
    return hostArrayOf(items.size(), [&](size_t i) { return items[i]; });
}

} // namespace bromesh::api

