#pragma once

// Adhesion model of Fang and Li et al. [2023].

#include <ipc/friction/smooth_mu.hpp>
#include <ipc/math/scalar_math.hpp>
#include <ipc/utils/simd.hpp>

#include <cassert>

namespace ipc {

// See the note atop smooth_friction_mollifier.hpp for why the branches here are
// written as `select_lazy` and what that means for a batch scalar.

// Fang and Li et al. [2023]:

// -- Normal Adhesion ----------------------------------------------------------
/// @defgroup normal_adhesion Normal Adhesion
/// @{

/// @brief The normal adhesion potential.
/// @tparam T The scalar type.
/// @param d distance
/// @param dhat_p distance of largest adhesion force (\f$\hat{d}_p\f$) where \f$0 < \hat{d}_p < \hat{d}_a\f$
/// @param dhat_a adhesion activation distance (\f$\hat{d}_a\f$)
/// @param a2 adjustable parameter relating to the maximum derivative of a (\f$a_2\f$)
/// @return The normal adhesion potential.
template <typename T>
inline T
normal_adhesion_potential(const T d, const T dhat_p, const T dhat_a, const T a2)
{
    assert(all_of(d >= T(0)));
    assert(all_of(dhat_p < dhat_a));
    assert(all_of(a2 < T(0)));
    return select_lazy(
        d < dhat_p,
        [&] {
            const T a1 = a2 * (T(1) - dhat_a / dhat_p);
            const T c1 = a2 * (dhat_a - dhat_p) * (dhat_a - dhat_p)
                - dhat_p * dhat_p * a1;
            return a1 * d * d + c1;
        },
        d < dhat_a,
        [&] {
            const T b2 = T(-2) * a2 * dhat_a;
            const T c2 = a2 * dhat_a * dhat_a;
            return (a2 * d + b2) * d + c2;
        },
        [&] { return T(0); });
}

/// @brief The first derivative of the normal adhesion potential wrt d.
/// @tparam T The scalar type.
/// @param d distance
/// @param dhat_p distance of largest adhesion force (\f$\hat{d}_p\f$) where \f$0 < \hat{d}_p < \hat{d}_a\f$
/// @param dhat_a adhesion activation distance (\f$\hat{d}_a\f$)
/// @param a2 adjustable parameter relating to the maximum derivative of a (\f$a_2\f$)
/// @return The first derivative of the normal adhesion potential wrt d.
template <typename T>
inline T normal_adhesion_potential_first_derivative(
    const T d, const T dhat_p, const T dhat_a, const T a2)
{
    assert(all_of(d >= T(0)));
    assert(all_of(dhat_p < dhat_a));
    assert(all_of(a2 < T(0)));
    return select_lazy(
        d < dhat_p,
        [&] {
            const T a1 = a2 * (T(1) - dhat_a / dhat_p);
            return T(2) * a1 * d;
        },
        d < dhat_a, [&] { return T(2) * a2 * (d - dhat_a); },
        [&] { return T(0); });
}

/// @brief The second derivative of the normal adhesion potential wrt d.
/// @tparam T The scalar type.
/// @param d distance
/// @param dhat_p distance of largest adhesion force (\f$\hat{d}_p\f$) where \f$0 < \hat{d}_p < \hat{d}_a\f$
/// @param dhat_a adhesion activation distance (\f$\hat{d}_a\f$)
/// @param a2 adjustable parameter relating to the maximum derivative of a (\f$a_2\f$)
/// @return The second derivative of the normal adhesion potential wrt d.
template <typename T>
inline T normal_adhesion_potential_second_derivative(
    const T d, const T dhat_p, const T dhat_a, const T a2)
{
    assert(all_of(d >= T(0)));
    assert(all_of(dhat_p < dhat_a));
    assert(all_of(a2 < T(0)));
    return select_lazy(
        d < dhat_p, [&] { return T(2) * a2 * (T(1) - dhat_a / dhat_p); },
        d < dhat_a, [&] { return T(2) * a2; }, [&] { return T(0); });
}

/// @brief The maximum normal adhesion force magnitude.
/// @tparam T The scalar type.
/// @param dhat_p distance of largest adhesion force (\f$\hat{d}_p\f$) where \f$0 < \hat{d}_p < \hat{d}_a\f$
/// @param dhat_a adhesion activation distance (\f$\hat{d}_a\f$)
/// @param a2 adjustable parameter relating to the maximum derivative of a (\f$a_2\f$)
/// @return The maximum normal adhesion force magnitude.
template <typename T>
inline T
max_normal_adhesion_force_magnitude(const T dhat_p, const T dhat_a, const T a2)
{
    assert(all_of(dhat_p < dhat_a));
    assert(all_of(a2 < T(0)));
    // max_d a' = a'(d̂ₚ) = 2a₂ (d̂ₚ - d̂ₐ)
    return T(2) * a2 * (dhat_p - dhat_a);
}

/// @}

// -- Tangential Adhesion ------------------------------------------------------
/// @defgroup tangential_adhesion Tangential Adhesion
/// @{

/// @brief The tangential adhesion mollifier function.
/// @tparam T The scalar type.
/// @param y The tangential relative speed.
/// @param eps_a Velocity threshold below which static adhesion force is applied.
/// @return The tangential adhesion mollifier function at y.
template <typename T> inline T tangential_adhesion_f0(const T y, const T eps_a)
{
    assert(all_of(eps_a > T(0)));
    return select_lazy(
        y <= T(0), [&] { return T(0); }, y >= T(2) * eps_a,
        [&] { return T(4) * eps_a / T(3); },
        // -y³/(3ϵ²) + y²/ϵ
        [&] { return y * y / eps_a * (T(1) - y / (T(3) * eps_a)); });
}

/// @brief The first derivative of the tangential adhesion mollifier function.
/// @tparam T The scalar type.
/// @param y The tangential relative speed.
/// @param eps_a Velocity threshold below which static adhesion force is applied.
/// @return The first derivative of the tangential adhesion mollifier function at y.
template <typename T> inline T tangential_adhesion_f1(const T y, const T eps_a)
{
    assert(all_of(eps_a > T(0)));
    return select_lazy(
        y >= T(2) * eps_a || y <= T(0), [&] { return T(0); },
        [&] {
            const T y_over_eps_a = y / eps_a;
            return y_over_eps_a * (T(2) - y_over_eps_a); // -y²/ϵ² + 2y/ϵ
        });
}

/// @brief The second derivative of the tangential adhesion mollifier function.
/// @tparam T The scalar type.
/// @param y The tangential relative speed.
/// @param eps_a Velocity threshold below which static adhesion force is applied.
/// @return The second derivative of the tangential adhesion mollifier function at y.
template <typename T> inline T tangential_adhesion_f2(const T y, const T eps_a)
{
    assert(all_of(eps_a > T(0)));
    return select_lazy(
        y >= T(2) * eps_a || y <= T(0), [&] { return T(0); },
        [&] { return (T(2) - T(2) * y / eps_a) / eps_a; }); // -2y/ϵ² + 2/ϵ
}

/// @brief The first derivative of the tangential adhesion mollifier function divided by y.
/// @tparam T The scalar type.
/// @param y The tangential relative speed.
/// @param eps_a Velocity threshold below which static adhesion force is applied.
/// @return The first derivative of the tangential adhesion mollifier function divided by y.
template <typename T>
inline T tangential_adhesion_f1_over_x(const T y, const T eps_a)
{
    assert(all_of(eps_a > T(0)));
    return select_lazy(
        y >= T(2) * eps_a || y <= T(0), [&] { return T(0); },
        [&] { return (T(2) - y / eps_a) / eps_a; }); // -y/ϵ² + 2/ϵ
}

/// @brief The second derivative of the tangential adhesion mollifier function times y minus the first derivative all divided by y cubed.
/// @tparam T The scalar type.
/// @param y The tangential relative speed.
/// @param eps_a Velocity threshold below which static adhesion force is applied.
/// @return The second derivative of the tangential adhesion mollifier function times y minus the first derivative all divided by y cubed.
template <typename T>
inline T tangential_adhesion_f2_x_minus_f1_over_x3(const T y, const T eps_a)
{
    assert(all_of(eps_a > T(0)));
    assert(all_of(y >= T(0)));
    return select_lazy(
        y >= T(2) * eps_a, [&] { return T(0); },
        [&] { return T(-1) / (y * eps_a * eps_a); });
}

// ~~ Smooth μ variants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// NOTE: Here a0, a1, and a2 refer to the mollifier functions above.

/// @brief Compute the value of the ∫ μ(y) a₁(y) dy, where a₁ is the first derivative of the smooth tangential adhesion mollifier.
/// @note The `a0`/`a1` are unrelated to the `a0`/`a1` in the normal adhesion.
/// @tparam T The scalar type.
/// @param y The tangential relative speed.
/// @param mu_s Coefficient of static adhesion.
/// @param mu_k Coefficient of kinetic adhesion.
/// @param eps_a Velocity threshold below which static adhesion force is applied.
/// @return The value of the integral at y.
template <typename T>
inline T smooth_mu_a0(const T y, const T mu_s, const T mu_k, const T eps_a)
{
    assert(all_of(eps_a > T(0)));
    const T delta_mu = mu_k - mu_s;
    const T z = y / eps_a;
    return select_lazy(
        y <= T(0), [&] { return T(0); }, //
        mu_s == mu_k || y >= eps_a,
        [&] {
            const T c = literal<T>(11 / 48.) * eps_a * delta_mu;
            return mu_k * tangential_adhesion_f0(y, eps_a) - c;
        },
        y < T(0.5) * eps_a,
        [&] {
            return y * z
                * (z
                       * (z * (T(1) - literal<T>(0.4) * z) * delta_mu
                          - mu_s / T(3))
                   + mu_s);
        },
        [&] {
            return y * z
                * (z
                       * (z * (literal<T>(0.4) * z - T(2)) * delta_mu
                          + T(3) * mu_k - literal<T>(10.0 / 3.0) * mu_s)
                   - mu_k + T(2) * mu_s)
                + literal<T>(3.0 / 80.0) * eps_a * delta_mu;
        });
}

/// @brief Compute the value of the μ(y) a₁(y), where a₁ is the first derivative of the smooth tangential adhesion mollifier.
/// @note The `a1` is unrelated to the `a1` in the normal adhesion.
/// @tparam T The scalar type.
/// @param y The tangential relative speed.
/// @param mu_s Coefficient of static adhesion.
/// @param mu_k Coefficient of kinetic adhesion.
/// @param eps_a Velocity threshold below which static adhesion force is applied.
/// @return The value of the product at y.
template <typename T>
inline T smooth_mu_a1(const T y, const T mu_s, const T mu_k, const T eps_a)
{
    return smooth_mu(y, mu_s, mu_k, eps_a) * tangential_adhesion_f1(y, eps_a);
}

/// @brief Compute the value of d/dy (μ(y) a₁(y)), where a₁ is the first derivative of the smooth tangential adhesion mollifier.
/// @note The `a1`/`a2` are unrelated to the `a1`/`a2` in the normal adhesion.
/// @tparam T The scalar type.
/// @param y The tangential relative speed.
/// @param mu_s Coefficient of static adhesion.
/// @param mu_k Coefficient of kinetic adhesion.
/// @param eps_a Velocity threshold below which static adhesion force is applied.
/// @return The value of the derivative at y.
template <typename T>
inline T smooth_mu_a2(const T y, const T mu_s, const T mu_k, const T eps_a)
{
    return smooth_mu_derivative(y, mu_s, mu_k, eps_a)
        * tangential_adhesion_f1(y, eps_a)
        + smooth_mu(y, mu_s, mu_k, eps_a) * tangential_adhesion_f2(y, eps_a);
}

/// @brief Compute the value of the μ(y) a₁(y) / y, where a₁ is the first derivative of the smooth tangential adhesion mollifier.
/// @note The `x` in the function name refers to the parameter `y`.
/// @note The `a1` is unrelated to the `a1` in the normal adhesion.
/// @tparam T The scalar type.
/// @param y The tangential relative speed.
/// @param mu_s Coefficient of static adhesion.
/// @param mu_k Coefficient of kinetic adhesion.
/// @param eps_a Velocity threshold below which static adhesion force is applied.
/// @return The value of the product at y.
template <typename T>
inline T
smooth_mu_a1_over_x(const T y, const T mu_s, const T mu_k, const T eps_a)
{
    // This is a known formulation: μ(y) f₁(y) / y
    // where we use the robust division by y to avoid division by zero.
    return smooth_mu(y, mu_s, mu_k, eps_a)
        * tangential_adhesion_f1_over_x(y, eps_a);
}

/// @brief Compute the value of the [(d/dy μ(y) a₁(y)) ⋅ y - μ(y) a₁(y)] / y³, where a₁ and a₂ are the first and second derivatives of the smooth tangential adhesion mollifier.
/// @note The `x` in the function name refers to the parameter `y`.
/// @note The `a1`/`a2` are unrelated to the `a1`/`a2` in the normal adhesion.
/// @tparam T The scalar type.
/// @param y The tangential relative speed.
/// @param mu_s Coefficient of static adhesion.
/// @param mu_k Coefficient of kinetic adhesion.
/// @param eps_a Velocity threshold below which static adhesion force is applied.
/// @return The value of the expression at y.
template <typename T>
inline T smooth_mu_a2_x_minus_mu_a1_over_x3(
    const T y, const T mu_s, const T mu_k, const T eps_a)
{
    assert(all_of(eps_a > T(0)));
    assert(all_of(y >= T(0)));
    const T delta_mu = mu_k - mu_s;
    const T z = T(1) / eps_a;
    return select_lazy(
        mu_s == mu_k || y >= eps_a,
        [&] {
            return mu_k * tangential_adhesion_f2_x_minus_f1_over_x3(y, eps_a);
        },
        y < T(0.5) * eps_a,
        [&] {
            return z * z * (z * (T(8) - T(6) * z * y) * delta_mu - mu_s / y);
        },
        [&] {
            return z * z
                * (z * (T(6) * z * y - T(16)) * delta_mu
                   + (T(9) * mu_k - T(10) * mu_s) / y);
        });
}

/// @}

} // namespace ipc
