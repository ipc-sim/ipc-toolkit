#pragma once

#include "smooth_friction_mollifier.hpp"

namespace ipc {

template <typename T> inline T f0_SF(const T& x, const double epsv_times_h)
{
    assert(epsv_times_h > 0);
    if (std::abs(x) >= epsv_times_h) {
        return x;
    }
    return x * x * (-x / (3 * epsv_times_h) + 1) / epsv_times_h
        + epsv_times_h / 3;
}

template <typename T>
inline T f1_SF_over_x(const T& x, const double epsv_times_h)
{
    assert(epsv_times_h > 0);
    if (std::abs(x) >= epsv_times_h) {
        return 1 / x;
    }
    return (-x / epsv_times_h + 2) / epsv_times_h;
}

template <typename T>
inline T df1_x_minus_f1_over_x3(const T& x, const double epsv_times_h)
{
    assert(epsv_times_h > 0);
    if (std::abs(x) >= epsv_times_h) {
        return -1 / (x * x * x);
    }
    return -1 / (x * epsv_times_h * epsv_times_h);
}

} // namespace ipc
