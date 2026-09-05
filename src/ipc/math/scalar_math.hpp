#pragma once

#include <ipc/config.hpp>

#include <cmath>
#include <type_traits>

// When compiling CUDA device code with NVCC pull in math functions from the
// global namespace. In host mode, and when device code is compiled with clang,
// use the std versions.
#if defined(IPC_TOOLKIT_WITH_CUDA) && defined(__NVCC__)
#define IPC_TOOLKIT_USING_STD(FUNC) using ::FUNC;
#else
#define IPC_TOOLKIT_USING_STD(FUNC) using std::FUNC;
#endif

namespace ipc {

/// @brief Define `ipc::FUNC` forwarding its arguments to the `std` counterpart.
#define IPC_TOOLKIT_DEFINE_STD(FUNC)                                           \
    template <                                                                 \
        typename T, typename... Ts,                                            \
        typename = std::enable_if_t<(std::is_same_v<T, Ts> && ...)>>           \
    inline auto FUNC(const T& x, const Ts&... rest)                            \
    {                                                                          \
        IPC_TOOLKIT_USING_STD(FUNC)                                            \
        return FUNC(x, rest...);                                               \
    }

IPC_TOOLKIT_DEFINE_STD(abs)
IPC_TOOLKIT_DEFINE_STD(atan2)
IPC_TOOLKIT_DEFINE_STD(fma)
IPC_TOOLKIT_DEFINE_STD(log)
IPC_TOOLKIT_DEFINE_STD(sqrt)

// This is a public header, so the macro does not outlive its use here.
#undef IPC_TOOLKIT_DEFINE_STD

constexpr double MOLLIFIER_THRESHOLD_EPS = 1e-2;

/// @brief Square of `x`, for any scalar the library templates on.
/// @note Faster than `std::pow(x, 2)`.
template <typename T> inline T sqr(const T& x) { return x * x; }

/// @brief Cube of `x`, for any scalar the library templates on.
template <typename T> inline T cubic(const T& x) { return x * x * x; }

} // namespace ipc
