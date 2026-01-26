// NOTE: This is an internal header file, not meant to be used outside of the
// IPC Toolkit library. It includes TBB which is a private dependency of the IPC
// Toolkit library. To use this outside of the library, one needs to link
// against TBB::tbb.

#pragma once

#include <ipc/utils/profiler.hpp>
#include <ipc/utils/unordered_map_and_set.hpp>

#include <tbb/enumerable_thread_specific.h>

#include <iterator> // std::make_move_iterator
#include <type_traits>
#include <vector>

namespace ipc {

// Assumes `out` is empty at the start. The function may modify the provided
// `vectors` (stealing and clearing per-thread buffers) for performance.
template <typename T>
void merge_thread_local_vectors(
    tbb::enumerable_thread_specific<std::vector<T>>& vectors,
    std::vector<T>& out)
{
    IPC_TOOLKIT_PROFILE_BLOCK("merge_thread_local_vectors");

    assert(out.empty());

    // Since `out` is always empty, compute total from thread-local vectors
    // only.
    size_t total = 0;
    for (auto& v : vectors) {
        total += v.size();
    }
    if (total == 0) {
        return;
    }

    // Fast path for trivially-copyable types: allocate once and memcpy each
    // thread-local buffer into the contiguous destination.
    if constexpr (
        std::is_trivially_copyable_v<T> && std::is_default_constructible_v<T>) {
        out.resize(total);
        char* dest = reinterpret_cast<char*>(out.data());
        for (auto& v : vectors) {
            if (v.empty()) {
                continue;
            }
            std::memcpy(dest, v.data(), v.size() * sizeof(T));
            dest += v.size() * sizeof(T);
            // release the local buffer to reduce memory usage
            std::vector<T>().swap(v);
        }
    } else {
        // For non-trivial types, steal the largest thread-local buffer into
        // `out` (cheap swap) and move from the remaining buffers.
        std::vector<T>* largest = nullptr;
        for (auto& v : vectors) {
            if (!largest || v.size() > largest->size()) {
                largest = &v;
            }
        }

        if (largest && !largest->empty()) {
            // out is empty, so swapping moves the largest contents into out and
            // leaves the thread-local buffer empty (former out).
            out.swap(*largest);
        }

        out.reserve(total);

        for (auto& v : vectors) {
            if (&v != largest && !v.empty()) {
                // Move elements to `out` to avoid copies when possible.
                out.insert(
                    out.end(), std::make_move_iterator(v.begin()),
                    std::make_move_iterator(v.end()));
            }
            // Ensure capacity is released.
            std::vector<T>().swap(v);
        }
    }
}

// No trait detection: always use the generic insert-based fallback for
// portability and simplicity.

template <typename T>
void merge_thread_local_unordered_sets(
    tbb::enumerable_thread_specific<unordered_set<T>>& sets,
    unordered_set<T>& out)
{
    IPC_TOOLKIT_PROFILE_BLOCK("merge_thread_local_unordered_sets");

    // This function assumes `out` is empty at the start and is allowed to
    // modify the per-thread `sets` (stealing / clearing buffers) for better
    // performance.
    assert(out.empty());

    // Compute total number of elements across thread-local sets.
    size_t total = 0;
    for (auto& s : sets) {
        total += s.size();
    }
    if (total == 0) {
        return;
    }

    // Steal the largest set by swapping it with `out` (cheap).
    unordered_set<T>* largest = nullptr;
    for (auto& s : sets) {
        if (!largest || s.size() > largest->size()) {
            largest = &s;
        }
    }

    if (largest && !largest->empty()) {
        out.swap(*largest);
    }

    out.reserve(total);

    // Simplified strategy: always insert remaining elements from per-thread
    // sets into `out`. After inserting we free the local buffer to reduce
    // memory usage. This avoids complex trait detection and keeps behavior
    // portable across different unordered_set implementations.
    for (auto& s : sets) {
        if (&s != largest && !s.empty()) {
            out.insert(s.begin(), s.end());
        }
        unordered_set<T>().swap(s);
    }
}

} // namespace ipc