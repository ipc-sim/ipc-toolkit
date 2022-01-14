#pragma once

#include <ipc/collision_constraint.hpp>
#include <ipc/friction/relative_displacement.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

// C1 clamping
template <typename T> inline T f0_SF(const T& x, const double& epsv_times_h)
{
    assert(epsv_times_h >= 0);
    if (abs(x) >= epsv_times_h) {
        return x;
    }
    return x * x * (-x / (3 * epsv_times_h) + 1) / epsv_times_h
        + epsv_times_h / 3;
}

/// Derivative of f0_SF divided by x
template <typename T>
inline T f1_SF_over_x(const T& x, const double& epsv_times_h)
{
    assert(epsv_times_h >= 0);
    if (abs(x) >= epsv_times_h) {
        return 1 / x;
    }
    return (-x / epsv_times_h + 2) / epsv_times_h;
}

template <typename T> inline T f2_SF(const T& x, const double& epsv_times_h)
{
    // same for abs(x) >= epsv_times_h for C1 clamped friction
    return T(-1 / (epsv_times_h * epsv_times_h));
}

} // namespace ipc
