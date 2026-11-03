// NOTE: This is an internal header file, not meant to be used outside of the
// IPC Toolkit library. It includes TBB which is a private dependency of the IPC
// Toolkit library. To use this outside of the library, one needs to link
// against TBB::tbb.

#pragma once

#include <ipc/utils/profiler.hpp>
#include <ipc/utils/unordered_map_and_set.hpp>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include <cstring>  // std::memcpy
#include <iterator> // std::make_move_iterator
#include <type_traits>
#include <vector>

namespace ipc {

/// @brief Payload size below which merging in parallel costs more than it
/// saves: the copies are short enough that the parallel_for's own overhead
/// dominates.
inline constexpr size_t PARALLEL_MERGE_MIN_BYTES = 1 << 20; // 1 MiB

// Assumes `out` is empty at the start. The function may modify the provided
// `vectors` (stealing and clearing per-thread buffers) for performance.
template <typename T, typename AllocIn, typename AllocOut>
void merge_thread_local_vectors(
    tbb::enumerable_thread_specific<std::vector<T, AllocIn>>& vectors,
    std::vector<T, AllocOut>& out)
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
    // thread-local buffer into the contiguous destination. Each buffer lands
    // in its own slice of `out`, so once the offsets are known with a prefix
    // sum the copies are independent and can run concurrently.
    if constexpr (
        std::is_trivially_copyable_v<T> && std::is_default_constructible_v<T>) {
        out.resize(total);

        if (total * sizeof(T) < PARALLEL_MERGE_MIN_BYTES) {
            // Copy in order. The bookkeeping the parallel path needs costs
            // two allocations, which is more than it saves at this size.
            char* dest = reinterpret_cast<char*>(out.data());
            for (auto& v : vectors) {
                if (v.empty()) {
                    continue;
                }
                std::memcpy(dest, v.data(), v.size() * sizeof(T));
                dest += v.size() * sizeof(T);
                // release the local buffer to reduce memory usage
                std::vector<T, AllocIn>().swap(v);
            }
        } else {
            std::vector<std::vector<T, AllocIn>*> blocks;
            blocks.reserve(vectors.size());
            for (auto& v : vectors) {
                if (!v.empty()) {
                    blocks.push_back(&v);
                }
            }

            std::vector<size_t> offsets(blocks.size() + 1, 0);
            for (size_t i = 0; i < blocks.size(); i++) {
                offsets[i + 1] = offsets[i] + blocks[i]->size();
            }

            tbb::parallel_for(size_t(0), blocks.size(), [&](const size_t i) {
                std::memcpy(
                    out.data() + offsets[i], blocks[i]->data(),
                    blocks[i]->size() * sizeof(T));
                // release the local buffer to reduce memory usage
                std::vector<T, AllocIn>().swap(*blocks[i]);
            });
        }
    } else {
        // For non-trivial types, steal the largest thread-local buffer into
        // `out` (cheap swap) and move from the remaining buffers.
        std::vector<T, AllocIn>* largest = nullptr;
        for (auto& v : vectors) {
            if (!largest || v.size() > largest->size()) {
                largest = &v;
            }
        }

        if constexpr (std::is_same_v<AllocIn, AllocOut>) {
            if (largest && !largest->empty()) {
                // out is empty, so swapping moves the largest contents into
                // out and leaves the thread-local buffer empty (former out).
                out.swap(*largest);
            }
        } else {
            // The buffers do not share an allocator, so the largest cannot be
            // stolen; it is moved from with the rest below.
            largest = nullptr;
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
            std::vector<T, AllocIn>().swap(v);
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