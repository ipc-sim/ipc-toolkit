#pragma once

#include <ipc/config.hpp>

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {
    template <typename scalar>
    scalar smooth_heaviside(const scalar &x);

    template <typename scalar>
    scalar heaviside(const scalar &x);

    template <typename scalar>
    scalar smooth_max(const scalar &x, const scalar &y);

    template <typename scalar>
    scalar smooth_min(const scalar &x, const scalar &y);

    template <typename scalar>
    scalar cubic_spline(const scalar &x);

    template <typename scalar>
    scalar inv_barrier(const scalar &x, const double eps, const double r);
}