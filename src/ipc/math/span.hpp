#include <cstddef> // for std::size_t

namespace ipc {
// A minimal, non-owning view of a contiguous sequence of objects.
// Named to mirror std::span, hence the non-CamelCase name.
// NOLINTNEXTLINE(readability-identifier-naming)
template <typename T> class span {
public:
    // Member types
    using element_type = T;
    using value_type = std::remove_cv_t<T>;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using pointer = T*;
    using const_pointer = const T*;
    using reference = T&;
    using const_reference = const T&;
    using iterator = pointer;
    using const_iterator = const_pointer;

    // Constructors
    // Default constructor (creates an empty span)
    constexpr span() noexcept : m_ptr(nullptr), m_size(0) { }

    // Construct from a pointer and a count
    constexpr span(pointer ptr, size_type count) noexcept
        : m_ptr(ptr)
        , m_size(count)
    {
    }

    // Construct from a pointer and an end pointer
    constexpr span(pointer first, pointer last) noexcept
        : m_ptr(first)
        , m_size(static_cast<size_type>(last - first))
    {
    }

    // Element access
    constexpr reference operator[](size_type idx) const noexcept
    {
        // In a real implementation, bounds checking might be optional (e.g., in
        // debug builds).
        return *(m_ptr + idx);
    }

    constexpr pointer data() const noexcept { return m_ptr; }

    // Observers
    constexpr size_type size() const noexcept { return m_size; }
    constexpr bool empty() const noexcept { return m_size == 0; }

    // Iterators
    constexpr iterator begin() const noexcept { return m_ptr; }
    constexpr iterator end() const noexcept { return m_ptr + m_size; }

private:
    pointer m_ptr;
    size_type m_size;
};
} // namespace ipc
