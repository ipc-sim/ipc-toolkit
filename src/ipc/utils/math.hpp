#pragma once

#include <ipc/config.hpp>

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {
    template <typename scalar>
    scalar intpow(const scalar &x, int p)
    {
        assert(p >= 0);
        if (p == 0)
            return scalar(1.);
        
        scalar out = x;
        while (p-- > 1)
            out = out * x;
        
        return out;
    }
    
    template <typename scalar>
    scalar cubic_spline(const scalar &x)
    {
        if (x <= -2)
            return scalar(0.);
        if (x <= -1)
            return intpow(x + 2, 3) / 6;
        if (x <= 0)
            return 2. / 3 - intpow(x, 2) * (1 + x / 2);
        if (x <= 1)
            return 2. / 3 - intpow(x, 2) * (1 - x / 2);
        if (x < 2)
            return -intpow(x - 2, 3) / 6;
        return scalar(0.);
    }

    template <typename scalar>
    scalar inv_barrier(const scalar &x, const double &eps, const double &r)
    {
        return cubic_spline((2 / eps) * x) / pow(x, r / 2.);
    }

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
    scalar L_s(const scalar &x, const double &a)
    {
        if (x < -a)
            return x;
        if (x < 0)
            return -intpow(x / a, 2) * (2 * a + x);
        if (x < 1)
            return scalar(0.);
        
        const scalar z = x - 1;
        if (z < a)
            return intpow(z / a, 2) * (2 * a - z);
        return z;
    }

    template <typename scalar>
    scalar cross2(const Eigen::Ref<const VectorMax3<scalar>> &a, const Eigen::Ref<const VectorMax3<scalar>> &b)
    {
        assert(a.size() == 3 || a.size() == 2);
        assert(b.size() == 3 || b.size() == 2);
        if (a.size() == 2)
            return a[0] * b[1] - a[1] * b[0];
        else if (a.size() == 3)
        {
            const Eigen::Ref<const Vector3<scalar>> a_(a), b_(b);
            return a_.cross(b_).norm();
        }
        else
            return scalar(0);
    }
}