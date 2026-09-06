#pragma once

#include <ipc/math/scalar_math.hpp>
#include <ipc/utils/simd.hpp>

#include <cassert>

namespace ipc {

// Each function below is piecewise around |y| = eps_v. We write the branch as
// `select_lazy` rather than an `if`, which is what lets one definition serve a
// plain `double`, an autodiff scalar, and a SIMD batch.
//
// The two paths differ in a way worth knowing before editing these: for a
// scalar, `select_lazy` is an ordinary `if`/`else` and only the winning branch
// runs. For a batch, lanes may straddle the threshold, so *both* branches run
// and the results are blended per-lane. Several of the inactive branches below
// divide by `y`, so on a lane where `y == 0` that branch produces an infinity —
// which is harmless, because the blend is a bitwise per-lane select that
// discards it rather than arithmetic that would propagate it into a NaN.

/// @brief Smooth friction mollifier function.
///
/// \f\[
///     f_0(y)= \begin{cases}
///         -\frac{y^3}{3\epsilon_v^2} + \frac{y^2}{\epsilon_v}
///             + \frac{\epsilon_v}{3}, & |y| < \epsilon_v
///             \newline
///         y, & |y| \geq \epsilon_v
///     \end{cases}
/// \f\]
///
/// @tparam T The scalar type.
/// @param y The tangential relative speed.
/// @param eps_v Velocity threshold below which static friction force is applied.
/// @return The value of the mollifier function at y.
template <typename T> inline T smooth_friction_f0(const T y, const T eps_v)
{
    assert(all_of(eps_v > T(0)));
    return select_lazy(
        ipc::numext::abs(y) >= eps_v, [&] { return y; },
        [&] {
            return y * y * (T(1) - y / (T(3) * eps_v)) / eps_v + eps_v / T(3);
        });
}

/// @brief The first derivative of the smooth friction mollifier.
///
/// \f\[
///     f_1(y) = f_0'(y) = \begin{cases}
///         -\frac{y^2}{\epsilon_v^2}+\frac{2 y}{\epsilon_v}, & |y| < \epsilon_v
///         \newline 1, & |y| \geq \epsilon_v
///     \end{cases}
/// \f\]
///
/// @tparam T The scalar type.
/// @param y The tangential relative speed.
/// @param eps_v Velocity threshold below which static friction force is applied.
/// @return The value of the derivative of the smooth friction mollifier at y.
template <typename T> inline T smooth_friction_f1(const T y, const T eps_v)
{
    assert(all_of(eps_v > T(0)));
    return select_lazy(
        ipc::numext::abs(y) >= eps_v, [&] { return T(1); },
        [&] {
            const T y_over_eps_v = y / eps_v;
            return y_over_eps_v * (T(2) - y_over_eps_v);
        });
}

/// @brief The second derivative of the smooth friction mollifier.
///
/// \f\[
///     f_2(y) = f_0''(y) = \begin{cases}
///         -\frac{2 y}{\epsilon_v^2}+\frac{2}{\epsilon_v}, & |y| < \epsilon_v
///         \newline 0, & |y| \geq \epsilon_v
///     \end{cases}
/// \f\]
///
/// @tparam T The scalar type.
/// @param y The tangential relative speed.
/// @param eps_v Velocity threshold below which static friction force is applied.
/// @return The value of the second derivative of the smooth friction mollifier at y.
template <typename T> inline T smooth_friction_f2(const T y, const T eps_v)
{
    assert(all_of(eps_v > T(0)));
    return select_lazy(
        ipc::numext::abs(y) >= eps_v, [&] { return T(0); },
        [&] { return (T(2) - T(2) * y / eps_v) / eps_v; });
}

/// @brief Compute the derivative of the smooth friction mollifier divided by y (\f$\frac{f_0'(y)}{y}\f$).
///
/// \f\[
///     \frac{f_1(y)}{y} = \begin{cases}
///         -\frac{y}{\epsilon_v^2}+\frac{2}{\epsilon_v}, & |y| < \epsilon_v
///         \newline \frac{1}{y}, & |y| \geq \epsilon_v
///     \end{cases}
/// \f\]
///
/// @note The `x` in the function name refers to the parameter `y`.
/// @tparam T The scalar type.
/// @param y The tangential relative speed.
/// @param eps_v Velocity threshold below which static friction force is applied.
/// @return The value of the derivative of smooth_friction_f0 divided by y.
template <typename T>
inline T smooth_friction_f1_over_x(const T y, const T eps_v)
{
    assert(all_of(eps_v > T(0)));
    return select_lazy(
        ipc::numext::abs(y) >= eps_v, [&] { return T(1) / y; },
        [&] { return (T(2) - y / eps_v) / eps_v; });
}

/// @brief The derivative of f1 times y minus f1 all divided by y cubed.
///
/// \f\[
///     \frac{f_1'(y) y - f_1(y)}{y^3} = \begin{cases}
///         -\frac{1}{y \epsilon_v^2}, & |y| < \epsilon_v \newline
///         -\frac{1}{y^3}, & |y| \geq \epsilon_v
///     \end{cases}
/// \f\]
///
/// @note The `x` in the function name refers to the parameter `y`.
/// @tparam T The scalar type.
/// @param y The tangential relative speed.
/// @param eps_v Velocity threshold below which static friction force is applied.
/// @return The derivative of f1 times y minus f1 all divided by y cubed.
template <typename T>
inline T smooth_friction_f2_x_minus_f1_over_x3(const T y, const T eps_v)
{
    assert(all_of(eps_v > T(0)));
    return select_lazy(
        ipc::numext::abs(y) >= eps_v, [&] { return T(-1) / (y * y * y); },
        [&] { return T(-1) / (y * eps_v * eps_v); });
}

} // namespace ipc
