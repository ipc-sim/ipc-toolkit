#pragma once

#include <tbb/enumerable_thread_specific.h>

#include <vector>

namespace ipc {

template <typename T>
void merge_thread_local_vectors(
    const tbb::enumerable_thread_specific<std::vector<T>>& vectors,
    std::vector<T>& out)
{
    // size up the items
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

template <typename T>
void merge_thread_local_unordered_sets(
    const tbb::enumerable_thread_specific<unordered_set<T>>& sets,
    unordered_set<T>& out)
{
    // size up the items
    size_t size = out.size();
    for (const auto& set : sets) {
        size += set.size();
    }
    // serial merge!
    out.reserve(size);
    for (const auto& set : sets) {
        out.insert(set.begin(), set.end());
    }
}

} // namespace ipc