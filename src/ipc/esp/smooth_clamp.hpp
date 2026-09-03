#pragma once

#include <type_traits>

namespace ipc {

/// Width of the cubic-blend region at each end of [0, 1] for smooth_clamp01.
constexpr double kSmoothClampEps = 0.1;

namespace detail {
    template <typename T> double smooth_clamp_scalar(const T& x)
    {
        if constexpr (std::is_same_v<T, double>)
            return x;
        else
            return x.val;
    }
} // namespace detail

/// C^1 smooth saturation onto [0, 1].
///   f(x) = 0                       for x <= 0
///   f(x) = -x^3/eps^2 + 2 x^2/eps  for 0 <= x <= eps      (cubic blend)
///   f(x) = x                       for eps <= x <= 1-eps
///   reflected cubic blend          for 1-eps <= x <= 1
///   f(x) = 1                       for x >= 1
/// Continuous and C^1 globally; monotone on (0, eps) since
///   f'(x) = (x/eps) * (4 - 3 x/eps) > 0 for x in (0, eps).
template <typename T> T smooth_clamp01(const T& x)
{
    constexpr double eps = kSmoothClampEps;
    const double xv = detail::smooth_clamp_scalar(x);
    if (xv <= 0.0)
        return T(0.0);
    if (xv >= 1.0)
        return T(1.0);
    if (xv < eps) {
        return -x * x * x / (eps * eps) + 2.0 * x * x / eps;
    }
    if (xv > 1.0 - eps) {
        const T s = 1.0 - x;
        return 1.0 + s * s * s / (eps * eps) - 2.0 * s * s / eps;
    }
    return x;
}

/// C^1 smooth projection of a barycentric pair (u, v) (with implicit
/// w = 1 - u - v) onto the 2-simplex { (a, b, c) : a, b, c >= 0, a+b+c = 1 }.
/// Independently smooth-clamps each component to [0, 1] then renormalizes.
///
/// Properties:
///   - Sum of returned components is 1 (out has u + v + w == 1).
///   - Identity on the interior hexagon [eps, 1-eps]^3 (since smooth_clamp01
///     is identity there and renormalization is by 1).
///   - C1 globally: each smooth_clamp01 is C1 and the denominator
///     u_s + v_s + w_s >= eps > 0 because the input satisfies u + v + w = 1
///     so at least one component is >= 1/3 > eps.
///   - Output (u, v) lies in the closed triangle.
template <typename T>
void smooth_clamp_simplex(const T& u, const T& v, T& u_out, T& v_out)
{
    const T w = 1.0 - u - v;
    const T u_s = smooth_clamp01(u);
    const T v_s = smooth_clamp01(v);
    const T w_s = smooth_clamp01(w);
    const T inv_sum = 1.0 / (u_s + v_s + w_s);
    u_out = u_s * inv_sum;
    v_out = v_s * inv_sum;
}

} // namespace ipc
