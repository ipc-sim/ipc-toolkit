#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_ABSEIL
#include <absl/hash/hash.h>
#else
#include <cstddef>
#include <functional>

namespace ipc {
template <typename H, typename T> H AbslHashValue(H h, T t);

template <class T> struct Hash {
    Hash() = default;
    Hash(size_t h) : hash(h) {};

    template <typename Value>
    static Hash combine(const Hash& h, const Value value)
    {
        if constexpr (std::is_default_constructible<std::hash<Value>>::value) {
            std::hash<Value> hash;
            return Hash(
                h.hash
                ^ (hash(value) + 0x9e3779b9 + (h.hash << 6) + (h.hash >> 2)));
        } else {
            return AbslHashValue(h, value);
        }
    }

    template <class First, class... Rest>
    static Hash combine(const Hash& h, const First first, const Rest... rest)
    {
        if constexpr (sizeof...(Rest) == 0) {
            return Hash::combine<First>(std::move(h), first);
        } else {
            return Hash::combine(Hash::combine(std::move(h), first), rest...);
        }
    }

    size_t operator()(const T& t) const
    {
        return AbslHashValue<Hash>(*this, t);
    }

    operator size_t() const { return hash; }

    size_t hash = 0;
};
} // namespace ipc
#endif

#ifdef IPC_TOOLKIT_WITH_ROBIN_MAP

#include <tsl/robin_map.h>
#include <tsl/robin_set.h>

namespace ipc {

#ifdef IPC_TOOLKIT_WITH_ABSEIL
template <typename K, typename V, typename Hash = absl::Hash<K>>
using unordered_map = tsl::robin_map<K, V, Hash>;
template <typename T, typename Hash = absl::Hash<T>>
using unordered_set = tsl::robin_set<T, Hash>;
#else
template <typename K, typename V, typename Hash = Hash<K>>
using unordered_map = tsl::robin_map<K, V, Hash>;
template <typename T, typename Hash = Hash<T>>
using unordered_set = tsl::robin_set<T, Hash>;
#endif

} // namespace ipc

#else

#include <unordered_map>
#include <unordered_set>

namespace ipc {

#ifdef IPC_TOOLKIT_WITH_ABSEIL
template <typename K, typename V, typename Hash = absl::Hash<K>>
using unordered_map = std::unordered_map<K, V, Hash>;
template <typename T, typename Hash = absl::Hash<T>>
using unordered_set = std::unordered_set<T, Hash>;
#else
template <typename K, typename V, typename Hash = Hash<K>>
using unordered_map = std::unordered_map<K, V, Hash>;
template <typename T, typename Hash = Hash<T>>
using unordered_set = std::unordered_set<T, Hash>;
#endif

} // namespace ipc

#endif
