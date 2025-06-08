#include "unordered_map_and_set.hpp"

#ifndef IPC_TOOLKIT_WITH_ABSEIL

#include <utility> // std::pair

namespace ipc {

template <>
Hash<std::pair<int, int>>
AbslHashValue(Hash<std::pair<int, int>> h, const std::pair<int, int> p)
{
    return Hash<std::pair<int, int>>::combine(std::move(h), p.first, p.second);
}

template <> Hash<int> AbslHashValue(Hash<int> h, const int i)
{
    return Hash<int>::combine(std::move(h), i);
}

template <>
Hash<std::pair<size_t, size_t>> AbslHashValue(
    Hash<std::pair<size_t, size_t>> h, const std::pair<size_t, size_t> p)
{
    return Hash<std::pair<size_t, size_t>>::combine(
        std::move(h), p.first, p.second);
}

template <> Hash<size_t> AbslHashValue(Hash<size_t> h, const size_t i)
{
    return Hash<size_t>::combine(std::move(h), i);
}

} // namespace ipc

#endif
