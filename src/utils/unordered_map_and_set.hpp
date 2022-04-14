#pragma once

#ifdef IPC_TOOLKIT_USE_ROBIN_MAP

#include <tsl/robin_map.h>
#include <tsl/robin_set.h>

namespace ipc {

template <typename K, typename V> using unordered_map = tsl::robin_map<K, V>;
template <typename T> using unordered_set = tsl::robin_set<T>;

} // namespace ipc

#else

#include <unordered_map>
#include <unordered_set>

namespace ipc {

using std::unordered_map;
using std::unordered_set;

} // namespace ipc

#endif
