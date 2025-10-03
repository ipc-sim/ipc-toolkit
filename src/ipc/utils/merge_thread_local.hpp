// NOTE: This is an internal header file, not meant to be used outside of the
// IPC Toolkit library. It includes TBB which is a private dependency of the IPC
// Toolkit library. To use this outside of the library, one needs to link
// against TBB::tbb.

#pragma once

#include <ipc/utils/unordered_map_and_set.hpp>

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