#pragma once

#ifdef IPC_TOOLKIT_WITH_TBB
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#else // Not using parallel for
#include <array>
#endif

namespace ipc {
namespace utils {

    template <typename T>
#ifdef IPC_TOOLKIT_WITH_TBB
    using ParallelCacheType = tbb::enumerable_thread_specific<T>;
#else
    using ParallelCacheType = std::array<T, 1>;
#endif

    // Perform a parallel (maybe) for loop.
    // The parallel for used depends on the compile definitions.
    // The overall for loop is from 0 up to `size` with an increment of 1.
    inline void maybe_parallel_for(
        int size, const std::function<void(int, int, int)>& partial_for);
    inline void
    maybe_parallel_for(int size, const std::function<void(int)>& body);

    // Returns thread specific storage for further use in
    // `maybe_parallel_for()`. The return type depends on the threading library
    // used.
    //     TBB         ⟹ `std::vector<LocalStorage>`
    //     C++ Threads ⟹ `tbb::enumerable_thread_specific<LocalStorage>`
    //     none        ⟹ `std::array<LocalStorage, 1>`
    template <typename LocalStorage>
    inline auto
    create_thread_storage(const LocalStorage& initial_local_storage);

    template <typename Storages>
    inline auto& get_local_thread_storage(Storages& storage, int thread_id);
} // namespace utils
} // namespace ipc

#include "maybe_parallel_for.tpp"
