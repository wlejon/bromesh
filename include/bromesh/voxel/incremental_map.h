#pragma once

#include <cstddef>
#include <unordered_map>
#include <utility>

namespace bromesh {

/// A hash map that grows without a pause: a std::unordered_map of half a
/// million entries rehashes them all in the insert that outgrows it, tens of
/// milliseconds in one call. Here that insert instead starts a table sized
/// for twice the entries, and every later insert or lookup moves a few
/// entries across from the old one (as node handles: no key is copied), so
/// the old table is empty well before the new one fills. Lookups look in
/// both while they coexist.
template <class K, class V, class H = std::hash<K>, class Eq = std::equal_to<K>>
class IncrementalMap {
public:
    std::size_t size() const noexcept { return cur_.size() + old_.size(); }

    V* find(const K& key) {
        step();
        if (auto it = cur_.find(key); it != cur_.end()) return &it->second;
        if (auto it = old_.find(key); it != old_.end()) return &it->second;
        return nullptr;
    }

    /// Insert a key that is not present.
    void emplace(K key, V value) {
        grow();
        step();
        cur_.emplace(std::move(key), std::move(value));
    }

    void erase(const K& key) {
        if (cur_.erase(key) == 0) old_.erase(key);
    }

    template <class Pred>
    void erase_if(Pred pred) {
        std::erase_if(cur_, pred);
        std::erase_if(old_, pred);
    }

    void clear() noexcept {
        cur_.clear();
        old_.clear();
    }

private:
    using Map = std::unordered_map<K, V, H, Eq>;
    // Below this a rehash is cheap: the table grows the ordinary way.
    static constexpr std::size_t kSmall = 4096;
    // Entries moved per call: the old table, n entries, is empty before the
    // new one (room for 2n) has taken n more.
    static constexpr int kMovesPerCall = 8;

    void grow() {
        if (cur_.size() < kSmall) return;
        if (static_cast<float>(cur_.size() + 1) <= cur_.max_load_factor() * static_cast<float>(cur_.bucket_count()))
            return;
        // Only reached with the old table drained (see kMovesPerCall); drain
        // what is left if not, so there are never three.
        while (!old_.empty()) cur_.insert(old_.extract(old_.begin()));
        old_ = std::move(cur_);
        cur_ = Map{};
        cur_.reserve(old_.size() * 2);
    }

    void step() {
        for (int i = 0; i < kMovesPerCall && !old_.empty(); ++i) cur_.insert(old_.extract(old_.begin()));
    }

    Map cur_;
    Map old_;
};

}  // namespace bromesh
