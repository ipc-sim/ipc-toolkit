#pragma once

namespace ipc {

/// @brief Square of `x`, for any scalar the library templates on.
/// @note Faster than `std::pow(x, 2)`.
template <typename T> inline T sqr(const T& x) { return x * x; }

/// @brief Cube of `x`, for any scalar the library templates on.
template <typename T> inline T cubic(const T& x) { return x * x * x; }

constexpr double MOLLIFIER_THRESHOLD_EPS = 1e-2;

} // namespace ipc
