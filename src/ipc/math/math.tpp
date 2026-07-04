#pragma once
#include "math.hpp"

#include <ipc/utils/autodiff_types.hpp>

namespace ipc {

namespace {
    /// @brief support is [0, 3]
    /// @tparam T
    /// @param x
    /// @return
    template <typename T> T quadratic_spline_aux(const T& x)
    {
        if (x <= 0) {
            return T(0.);
        }
        if (x <= 1) {
            return Math<T>::sqr(x) / 2.;
        }
        if (x <= 2) {
            return (-3 + x * (6 - 2 * x)) / 2.;
        }
        if (x < 3) {
            return Math<T>::sqr(3. - x) / 2.;
        }
        return T(0.);
    }

    template <typename T> T smooth_heaviside_standard(const T& x)
    {
        if (x <= -3) {
            return T(0.);
        }
        if (x <= -2) {
            return Math<T>::cubic(3. + x) / 6.;
        }
        if (x <= -1) {
            return (((-2 * x - 9) * x - 9) * x + 3) / 6;
        }
        if (x < 0) {
            return Math<T>::cubic(x) / 6. + 1.;
        }

        return T(1.);
    }

    [[maybe_unused]] double smooth_heaviside_standard_grad(const double x)
    {
        if (x <= -3 || x >= 0) {
            return 0.;
        }
        if (x <= -2) {
            return Math<double>::sqr(3. + x) / 2.;
        }
        if (x <= -1) {
            return -(x * x + 3 * x + 1.5);
        }

        return Math<double>::sqr(x) / 2.;
    }

    [[maybe_unused]] double smooth_heaviside_standard_hess(const double x)
    {
        if (x <= -3 || x >= 0) {
            return 0.;
        }
        if (x <= -2) {
            return 3. + x;
        }
        if (x <= -1) {
            return -3 - 2 * x;
        }

        return x;
    }
} // namespace

template <typename T> T Math<T>::cubic_spline(const T& x)
{
    if (abs(x) >= 1) {
        return T(0.);
    }
    if (abs(x) >= 0.5) {
        return cubic(1 - abs(x)) * (4. / 3.);
    }

    return 2. / 3. - 4. * (x * x) * (1 - abs(x));
}
template <typename T> double Math<T>::cubic_spline_grad(const double x)
{
    if (Math<double>::abs(x) >= 1) {
        return 0.;
    }
    if (Math<double>::abs(x) >= 0.5) {
        return -4. * Math<double>::sqr(1 - Math<double>::abs(x)) * sign(x);
    }

    return 4. * x * (3. * Math<double>::abs(x) - 2.);
}
template <typename T> double Math<T>::cubic_spline_hess(const double x)
{
    if (Math<double>::abs(x) >= 1) {
        return 0.;
    }
    if (Math<double>::abs(x) >= 0.5) {
        return 8. * (1 - Math<double>::abs(x));
    }

    return 8. * (3. * Math<double>::abs(x) - 1.);
}

/// @brief support is [-1, 1]
/// @tparam T
/// @param x
/// @return
template <typename T> T Math<T>::quadratic_spline(const T& x)
{
    return quadratic_spline_aux(x * 1.5 + 1.5);
}

// template <typename T>
// T smooth_heaviside_aux(const T &x)
// {
//     if (x <= -2)
//         return T(0.);
//     if (x <= -1)
//         return intpow(2. + x, 2) / 2.;
//     if (x <= 0)
//         return 1 - intpow(x, 2) / 2.;

//     return T(1.);
// }

template <typename T>
T Math<T>::smooth_heaviside(const T& x, const double alpha, const double beta)
{
    return smooth_heaviside_standard((x - beta) * (3 / (alpha + beta)));
}

template <typename T>
double Math<T>::smooth_heaviside_grad(
    const double x, const double alpha, const double beta)
{
    const double s = 3 / (alpha + beta);
    return smooth_heaviside_standard_grad((x - beta) * s) * s;
}

template <typename T>
double Math<T>::smooth_heaviside_hess(
    const double x, const double alpha, const double beta)
{
    const double s = 3 / (alpha + beta);
    return smooth_heaviside_standard_hess((x - beta) * s) * s * s;
}

template <typename T> T Math<T>::mollifier(const T& x)
{
    if constexpr (IsADHessian<T>::value) {
        if (x <= 0) {
            return T(0.);
        } else if (x < 1) {
            const double deriv = 2. * (1. - x.val), hess = -2.;
            return T::known_derivatives(
                x.val * (2. - x.val), deriv * x.grad,
                x.grad * hess * x.grad.transpose() + deriv * x.Hess);
        } else {
            return T(1.);
        }
        // return smooth_heaviside<T>(x - 1.);
    } else {
        if (x <= 0) {
            return T(0.);
        } else if (x < 1) {
            return x * (2. - x);
        } else {
            return T(1.);
        }
        // return smooth_heaviside<T>(x - 1.);
    }
}

template <typename T> double Math<T>::mollifier_grad(const double x)
{
    if (x <= 0 || x >= 1) {
        return 0.;
    } else {
        return 2. * (1. - x);
    }
}

template <typename T> double Math<T>::mollifier_hess(const double x)
{
    if (x <= 0 || x >= 1) {
        return 0.;
    } else {
        return -2.;
    }
}

// support is [0, 1]
template <typename T> T Math<T>::inv_barrier(const T& x, const int r)
{
    return cubic_spline(x) / pow(x, r);
    // log barrier
    // if (x < 1)
    //     return -(1 - x) * (1 - x) * log(x);
    // else
    //     return T(0.);
}

template <typename T>
double Math<T>::inv_barrier_grad(const double x, const int r)
{
    return (cubic_spline_grad(x) - Math<double>::cubic_spline(x) * r / x)
        / pow(x, r);
}
template <typename T>
double Math<T>::inv_barrier_hess(const double x, const int r)
{
    return (cubic_spline_hess(x)
            + (-2. * cubic_spline_grad(x)
               + (r + 1.) * Math<double>::cubic_spline(x) / x)
                * r / x)
        / pow(x, r);
}

template <typename T> T Math<T>::l_ns(const T& x)
{
    if (x <= 0.) {
        return T(0.);
    }
    if (x >= 1.) {
        return T(1.);
    }
    return x;
}

} // namespace ipc