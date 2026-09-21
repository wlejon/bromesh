#pragma once

#include "embed/embed.h"
#include "object_builder.h"

#include <functional>
#include <memory>
#include <span>
#include <string>

namespace bromesh::api {

namespace ev = bronze::embed;
using Value = bronze::Value;

// A host class: one constructor function + one prototype, with instances
// born on that prototype by make()/createInstance().
//
// The class objects themselves are process-global (`HostClass g_meshClass`),
// but what they hold is PER THREAD: bronze's runtime is per-thread, and a
// Persistent is a slot in its creating thread's registry, so a constructor
// made on the main thread means nothing to a Worker's realm. Every accessor
// below therefore reads the CALLING thread's slots, and install() fills the
// calling thread's. A realm installs each class once — the main realm from
// the host's sibling installer, a Worker realm from its own — and
// installed() answers for the calling thread, which is what an install
// guard has to ask.
class HostClass {
public:
    struct Slots {
        ev::Persistent* proto = nullptr;
        ev::Persistent* ctor = nullptr;
        ev::Persistent* instanceProto = nullptr;
        std::vector<std::string> aliases;
    };

    void install(const char* name, uint32_t arity, ev::NativeFn body,
                 const std::function<void(ObjectBuilder&)>& decorate = nullptr);

    void init(const char* name, const std::function<void(ObjectBuilder&)>& decorate) {
        install(name, 0, nullptr, decorate);
    }

    template <typename T>
    Value createInstance(std::unique_ptr<T> ptr) const {
        return make(ptr.release(), [](void* p) { delete static_cast<T*>(p); });
    }

    void alias(const char* name) const;
    void inherit(const HostClass& base) const;

    Value make(void* data, ev::HandleDestructor dtor,
               ev::Finalize when = ev::Finalize::InSweep) const;

    void setStatic(const char* name, Value v) const;

    void* unwrap(Value val) const { return ev::handleData(val); }

    Value prototype() const;
    Value constructor() const;
    void setInstancePrototype(Value p);

    // Whether install() has run on the CALLING thread.
    bool installed() const;

private:
    // The calling thread's slots: created empty on first touch.
    Slots& slots() const;
    // The calling thread's slots, or null when this thread never installed.
    const Slots* slotsIfAny() const;
};

} // namespace bromesh::api
