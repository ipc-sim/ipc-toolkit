#pragma once

#include <memory>

namespace ipc {

// A wrapper that performs no initialization on construction
template <typename T, typename A = std::allocator<T>>
class DefaultInitAllocator : public A {
    typedef std::allocator_traits<A> a_t;

public:
    // NOLINTNEXTLINE(readability-identifier-naming)
    template <typename U> struct rebind {
        using other =
            DefaultInitAllocator<U, typename a_t::template rebind_alloc<U>>;
    };
    using A::A;

    // The magic happens here: do nothing when constructing a U
    template <typename U>
    void
    construct(U* ptr) noexcept(std::is_nothrow_default_constructible<U>::value)
    {
        ::new (static_cast<void*>(ptr)) U; // Default construction
    }

    // Override construct to do nothing for trivial types
    template <typename U> void construct(U* ptr, const U& val)
    {
        a_t::construct(static_cast<A&>(*this), ptr, val);
    }
    template <typename U, typename... Args>
    void construct(U* ptr, Args&&... args)
    {
        a_t::construct(
            static_cast<A&>(*this), ptr, std::forward<Args>(args)...);
    }
};

} // namespace ipc