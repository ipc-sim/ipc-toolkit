#pragma once

#include <ipc/utils/default_init_allocator.hpp>

#include <vector>

namespace ipc {

/// @brief Storage for broad-phase candidates.
///
/// A std::vector that does not value-initialize on resize. The broad phase
/// sizes the output and then immediately overwrites every element, so zeroing
/// it first is a wasted pass over tens of megabytes: on Cloth-Ball's 5.2M
/// edge-edge candidates it costs 3.3 ms of a 3.8 ms merge.
///
/// The elements are left uninitialized by a resize that grows the vector, so
/// read only what has been written.
///
/// @tparam T The candidate type.
template <typename T>
using CandidateVector = std::vector<T, DefaultInitAllocator<T>>;

} // namespace ipc
