#pragma once

#include <ipc/config.hpp>

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {
    // template <typename scalar>
    // scalar smooth_heaviside(const scalar &x);

    // template <typename scalar>
    // scalar heaviside(const scalar &x);

    // template <typename scalar>
    // scalar smooth_max(const scalar &x, const scalar &y);

    // template <typename scalar>
    // scalar smooth_min(const scalar &x, const scalar &y);

    template <typename scalar>
    scalar cubic_spline(const scalar &x);

    template <typename scalar>
    scalar inv_barrier(const scalar &x, const double eps, const double r);

    template <typename scalar>
    scalar L_ns(const scalar &x)
    {
        if (x < 0)
            return x;
        if (x > 1)
            return x - 1;
        return scalar(0.);
    }

    template <typename scalar>
    scalar L_s(const scalar &x, int power)
    {
        return pow(L_ns(x), power);
    }

    template <typename scalar>
    scalar pow(const scalar &x, int p)
    {
        assert(p >= 1);
        scalar out = x;
        while (p-- > 1)
            out = out * x;
        
        return out;
    }
}