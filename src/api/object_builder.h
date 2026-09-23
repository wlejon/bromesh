#pragma once

#include "embed/embed.h"

#include <functional>
#include <span>
#include <string>
#include <string_view>
#include <utility>

namespace bromesh::api {

namespace ev = bronze::embed;
using Value = bronze::Value;

// ── Oversized lists ────────────────────────────────────────────────────
//
// The list readers that answer with a plain vector (toFloatVector and its
// siblings, host_mesh_internal.h) have no error channel, and a throw from
// inside one would not survive: the caller's next getProperty clears the
// pending exception at the host boundary. So a reader that meets a list
// longer than its cap (an array-like whose `length` claims billions of
// elements) records the refusal here and answers empty, and the binding
// boundary, hostFunction below, raises it as a RangeError once the binding
// returns, whatever the binding made of the empty list meanwhile.
inline thread_local std::string t_listRefusal;

inline void refuseList(std::string message) {
    if (t_listRefusal.empty()) t_listRefusal = std::move(message);
}

/// ev::makeFunction for every bromesh binding: the same function, with an
/// oversized list met anywhere inside it raised as a RangeError. A binding
/// that re-enters another through user code keeps its own record.
inline Value hostFunction(ev::NativeFn fn, uint32_t arity, std::string_view name) {
    ev::NativeFn guarded = [fn = std::move(fn)](Value self, std::span<const Value> args) -> Value {
        std::string outer = std::move(t_listRefusal);
        t_listRefusal.clear();
        Value r = fn(self, args);
        std::string mine = std::move(t_listRefusal);
        t_listRefusal = std::move(outer);
        if (!mine.empty()) return ev::throwRangeError(mine);
        return r;
    };
    return ev::makeFunction(std::move(guarded), arity, name);
}

/// Helper to build objects and namespaces property by property using bronze::embed.
/// Handles moving GC by rooting the target in an ev::Persistent.
struct ObjectBuilder {
    ev::Persistent obj;

    ObjectBuilder() : obj(ev::createObject()) {}
    explicit ObjectBuilder(Value existing) : obj(existing) {}

    void set(std::string_view name, Value v) {
        obj.set(ev::setProperty(obj.get(), name, v));
    }

    void set(std::string_view name, double d) {
        set(name, ev::fromDouble(d));
    }

    void set(std::string_view name, bool b) {
        set(name, ev::fromBool(b));
    }

    void set(std::string_view name, const std::string& s) {
        set(name, ev::fromUtf8(s));
    }

    void set(std::string_view name, const char* s) {
        set(name, ev::fromUtf8(s));
    }

    void def(std::string_view name, uint32_t arity, ev::NativeFn fn) {
        Value f = hostFunction(std::move(fn), arity, name);
        obj.set(ev::setProperty(obj.get(), name, f));
    }

    void accessor(std::string_view name, ev::NativeFn getter, ev::NativeFn setter = nullptr) {
        const std::string getName = "get " + std::string(name);
        const std::string setName = "set " + std::string(name);
        ev::Persistent g(hostFunction(std::move(getter), 0, getName));
        Value s = setter ? hostFunction(std::move(setter), 1, setName)
                         : ev::undefined();
        obj.set(ev::defineAccessor(obj.get(), name, g.get(), s, /*enumerable=*/true));
    }

    Value get() const { return obj.get(); }
    Value build() const { return obj.get(); }
};

} // namespace bromesh::api
