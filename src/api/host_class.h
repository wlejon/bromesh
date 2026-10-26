#pragma once

#include "embed/embed.h"
#include "embed/embed_host_class.h"
#include "object_builder.h"

#include <functional>
#include <memory>
#include <span>
#include <string>

namespace bromesh::api {

namespace ev = bronze::embed;
using Value = bronze::Value;

// bronze's HostClass (embed/embed_host_class.h): a constructor and prototype
// per realm, branded handles, and unwrap() checked across inherit() lineages.
// It used to be a copy kept per library, one constructor and prototype per
// thread; the realms of one thread then shared a class's mutable prototype,
// so one realm could change what another's instances do.
//
// The one thing added here is the decorate callback's builder: bromesh's own
// ObjectBuilder (object_builder.h), whose every method goes through
// hostFunction so an oversized list met inside it is raised as a RangeError.
// It decorates the prototype embed's install minted, and hands the object's
// final address back.
class HostClass : public bronze::embed::HostClass {
public:
    void install(const char* name, uint32_t arity, ev::NativeFn body,
                 const std::function<void(ObjectBuilder&)>& decorate = nullptr) {
        bronze::embed::HostClass::install(
            name, arity, std::move(body), [&decorate](ev::ObjectBuilder& proto) {
                if (!decorate) return;
                ObjectBuilder mine(proto.get());
                decorate(mine);
                proto.obj.set(mine.get());
            });
    }

    void init(const char* name, const std::function<void(ObjectBuilder&)>& decorate) {
        install(name, 0, nullptr, decorate);
    }
};

} // namespace bromesh::api
