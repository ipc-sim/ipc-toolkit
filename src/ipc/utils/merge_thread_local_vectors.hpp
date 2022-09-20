#pragma once

#include <tbb/enumerable_thread_specific.h>
#include <vector>

namespace ipc {

template <typename T>
void merge_thread_local_vectors(
    const tbb::enumerable_thread_specific<std::vector<T>>& vectors,
    std::vector<T>& out)
{
    // size up the hash items
    size_t size = out.size();
    for (const auto& vector : vectors) {
        size += vector.size();
    }
    // serial merge!
    out.reserve(size);
    for (const auto& vector : vectors) {
        out.insert(out.end(), vector.begin(), vector.end());
    }
}

} // namespace ipc