#pragma once

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

namespace bromesh {

/// An array that grows a page at a time and never moves what it holds: a
/// std::vector of millions of elements copies them all when it outgrows its
/// capacity, which for a voxel pool is tens of milliseconds in one call. Here
/// growing allocates one more page of `kPage` elements (each page contiguous,
/// so a GPU mirror uploads it page by page). Elements are value-initialised.
template <class T, std::size_t kPage = 4096>
class PagedArray {
public:
    static constexpr std::size_t kPageElems = kPage;

    std::size_t size() const noexcept { return size_; }
    bool empty() const noexcept { return size_ == 0; }

    T& operator[](std::size_t i) noexcept { return pages_[i / kPage][i % kPage]; }
    const T& operator[](std::size_t i) const noexcept { return pages_[i / kPage][i % kPage]; }
    T& back() noexcept { return (*this)[size_ - 1]; }
    const T& back() const noexcept { return (*this)[size_ - 1]; }

    template <class... A>
    T& emplace_back(A&&... args) {
        reserve(size_ + 1);
        T& slot = (*this)[size_];
        slot = T{std::forward<A>(args)...};
        ++size_;
        return slot;
    }
    void push_back(const T& v) { emplace_back(v); }

    /// Grows filled with `fill`, or shrinks (the pages stay).
    void resize(std::size_t n, const T& fill = T{}) {
        reserve(n);
        for (std::size_t i = size_; i < n; ++i) (*this)[i] = fill;
        size_ = n;
    }
    void clear() noexcept {
        pages_.clear();
        size_ = 0;
    }
    /// Pages for at least `n` elements.
    void reserve(std::size_t n) {
        while (pages_.size() * kPage < n) pages_.push_back(std::make_unique<T[]>(kPage));
    }

    /// The pages holding the elements, each `kPage` long but the last
    /// (`pageSize`).
    std::size_t pageCount() const noexcept { return (size_ + kPage - 1) / kPage; }
    const T* page(std::size_t p) const noexcept { return pages_[p].get(); }
    std::size_t pageSize(std::size_t p) const noexcept {
        const std::size_t first = p * kPage;
        return size_ - first < kPage ? size_ - first : kPage;
    }

private:
    std::vector<std::unique_ptr<T[]>> pages_;
    std::size_t size_ = 0;
};

}  // namespace bromesh
