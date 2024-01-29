#include "unordered_map_and_set.hpp"

#ifndef IPC_TOOLKIT_WITH_ABSEIL

#include <utility> // std::pair

namespace ipc {

template <>
Hash<std::pair<int, int>>
AbslHashValue(Hash<std::pair<int, int>> h, std::pair<int, int> p)
{
    return Hash<std::pair<int, int>>::combine(std::move(h), p.first, p.second);
}

template <>
Hash<std::pair<long, long>>
AbslHashValue(Hash<std::pair<long, long>> h, std::pair<long, long> p)
{
    return Hash<std::pair<long, long>>::combine(std::move(h), p.first, p.second);
}

template <> Hash<int> AbslHashValue(Hash<int> h, const int i)
{
    return Hash<int>::combine(std::move(h), i);
}

} // namespace ipc

#endif
