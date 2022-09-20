#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_USE_ROBIN_MAP

#include <tsl/robin_map.h>
#include <tsl/robin_set.h>

#ifdef IPC_TOOLKIT_USE_ABSL_HASH
#include <absl/hash/hash.h>
#endif

namespace ipc {

#ifdef IPC_TOOLKIT_USE_ABSL_HASH
template <typename K, typename V, typename Hash = absl::Hash<K>>
using unordered_map = tsl::robin_map<K, V, Hash>;
template <typename T, typename Hash = absl::Hash<T>>
using unordered_set = tsl::robin_set<T, Hash>;
#else
template <typename K, typename V> using unordered_map = tsl::robin_map<K, V>;
template <typename T> using unordered_set = tsl::robin_set<T>;
#endif

} // namespace ipc

#else

#include <unordered_map>
#include <unordered_set>

namespace ipc {

#ifdef IPC_TOOLKIT_USE_ABSL_HASH
template <typename K, typename V>
using unordered_map = std::unordered_map<K, V, absl::Hash<K>>;
template <typename T>
using unordered_set = std::unordered_set<T, absl::Hash<T>>;
#else
using std::unordered_map;
using std::unordered_set;
#endif

} // namespace ipc

#endif
