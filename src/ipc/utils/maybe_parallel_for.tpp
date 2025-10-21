#include "maybe_parallel_for.hpp"

#ifdef IPC_TOOLKIT_WITH_TBB
#include <tbb/parallel_for.h>
#endif

namespace ipc {

inline void maybe_parallel_for(
    int size, const std::function<void(int, int, int)>& partial_for)
{
#ifdef IPC_TOOLKIT_WITH_TBB
    tbb::parallel_for(
        tbb::blocked_range<int>(0, size),
        [&](const tbb::blocked_range<int>& r) {
            partial_for(
                r.begin(), r.end(),
                tbb::this_task_arena::current_thread_index());
        });
#else
    partial_for(0, size, /*thread_id=*/0); // actually the full for loop
#endif
}

inline void maybe_parallel_for(int size, const std::function<void(int)>& body)
{
#ifdef IPC_TOOLKIT_WITH_TBB
    tbb::parallel_for(0, size, body);
#else
    for (int i = 0; i < size; ++i) {
        body(i);
    }
#endif
}

template <typename LocalStorage>
inline auto create_thread_storage(const LocalStorage& initial_local_storage)
{
#ifdef IPC_TOOLKIT_WITH_TBB
    return tbb::enumerable_thread_specific<LocalStorage>(initial_local_storage);
#else
    return std::array<LocalStorage, 1> { { initial_local_storage } };
#endif
}

template <typename Storages>
inline auto& get_local_thread_storage(Storages& storage, int thread_id)
{
#ifdef IPC_TOOLKIT_WITH_TBB
    return storage.local();
#else
    assert(thread_id == 0);
    assert(storage.size() == 1);
    return storage[0];
#endif
}

} // namespace ipc
