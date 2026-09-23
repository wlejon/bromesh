#pragma once

#include "embed/embed.h"

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

inline int32_t i32At(std::span<const Value> args, size_t i) {
    return static_cast<int32_t>(static_cast<int64_t>(numAt(args, i)));
}

inline uint32_t u32At(std::span<const Value> args, size_t i) {
    return static_cast<uint32_t>(static_cast<int64_t>(numAt(args, i)));
}

inline int64_t i64At(std::span<const Value> args, size_t i) {
    return static_cast<int64_t>(numAt(args, i));
}

inline uint64_t u64At(std::span<const Value> args, size_t i) {
    if (i >= args.size()) return 0;
    Value v = args[i];
    if (ev::isObject(v)) return 0;
    return ev::toUint64(v);
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

