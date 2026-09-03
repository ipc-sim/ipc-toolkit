#pragma once

#include <cmath>

namespace ipc {

/// @brief Square of `x`, for any scalar the library templates on.
/// @note Faster than `std::pow(x, 2)`.
template <typename T> inline T sqr(const T& x) { return x * x; }

/// @brief Cube of `x`, for any scalar the library templates on.
template <typename T> inline T cubic(const T& x) { return x * x * x; }

/// @brief `sqrt` for any scalar the library templates on.
///
/// The block-scope using-declaration is what makes this work: unqualified
/// lookup stops at it (so this does not recurse), while ADL still reaches
/// `xsimd::sqrt` for a batch and TinyAD's hidden friend for an autodiff
/// scalar.
template <typename T> inline T sqrt(const T& x)
{
    using std::sqrt;
    return sqrt(x);
}

/// @brief `log` for any scalar the library templates on.
///
/// Same block-scope using-declaration trick as `sqrt` above: without it,
/// unqualified `log(float_val)` finds only the global `::log(double)`
/// declared by `<cmath>`, promoting the argument to `double` and narrowing
/// the result back, which trips `-Wfloat-conversion` for no reason.
template <typename T> inline T log(const T& x)
{
    using std::log;
    return log(x);
}

/// @brief `fma` for any scalar the library templates on.
///
/// Same block-scope using-declaration trick as `sqrt`/`log` above, but the
/// stakes are different: `std::fma` has no generic/arithmetic-type overload
/// that would silently compile for a batch, so a qualified `std::fma` call
/// on `SimdBatch<T>` is a hard compile error, not a quiet precision bug.
/// Unqualified lookup here lets ADL reach `xsimd::fma`.
template <typename T> inline T fma(const T& x, const T& y, const T& z)
{
    using std::fma;
    return fma(x, y, z);
}

constexpr double MOLLIFIER_THRESHOLD_EPS = 1e-2;

} // namespace ipc
