#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_TBB
#include <tbb/enumerable_thread_specific.h>
#else // Not using parallel for
#include <array>
#endif

/// @file maybe_parallel_for.hpp
/// @brief Utilities for conditional parallel for loops and thread-local storage.

namespace ipc {

/// @typedef ParallelCacheType
/// @brief Type alias for thread-local storage depending on threading library.
///
/// If IPC_TOOLKIT_WITH_TBB is defined, uses tbb::enumerable_thread_specific<T>.
/// Otherwise, uses std::array<T, 1>.
///
/// @tparam T Type of the thread-local storage.
template <typename T>
#ifdef IPC_TOOLKIT_WITH_TBB
using ParallelCacheType = tbb::enumerable_thread_specific<T>;
#else
using ParallelCacheType = std::array<T, 1>;
#endif

/// @brief Executes a (possibly parallel) for loop over a range.
///
/// The loop runs from 0 to `size` (exclusive) with an increment of 1.
/// The provided function is called with the start, end, and thread id for each
/// partition. The parallelization depends on compile-time definitions.
///
/// @param size Number of iterations.
/// @param partial_for Function to execute for each partition: void(int start, int end, int thread_id).
inline void maybe_parallel_for(
    int size, const std::function<void(int, int, int)>& partial_for);

/// @brief Executes a (possibly parallel) for loop over a range.
///
/// The loop runs from 0 to `size` (exclusive) with an increment of 1.
/// The provided function is called for each index.
/// The parallelization depends on compile-time definitions.
///
/// @param size Number of iterations.
/// @param body Function to execute for each index: void(int).
inline void maybe_parallel_for(int size, const std::function<void(int)>& body);

/// @brief Creates thread-specific storage for use in maybe_parallel_for().
///
/// The return type depends on the threading library:
/// - TBB: tbb::enumerable_thread_specific<LocalStorage>
/// - None: std::array<LocalStorage, 1>
///
/// @tparam LocalStorage Type of the local storage.
/// @param initial_local_storage Initial value for each thread's local storage.
/// @return Thread-specific storage container.
template <typename LocalStorage>
inline auto create_thread_storage(const LocalStorage& initial_local_storage);

/// @brief Retrieves the local storage for a specific thread.
///
/// @tparam Storages Type of the storage container.
/// @param storage Thread-specific storage container.
/// @param thread_id Identifier of the thread.
/// @return Reference to the local storage for the given thread.
template <typename Storages>
inline auto& get_local_thread_storage(Storages& storage, int thread_id);

} // namespace ipc

#include "maybe_parallel_for.tpp"
