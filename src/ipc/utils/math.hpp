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
    scalar mollifier(const scalar &x)
    {
        if (x <= 1)
            return x * x * (2. - x);
        return scalar(1.);
    }

    enum class HEAVISIDE_TYPE {
        ZERO = 0, ONE = 1, VARIANT = 2
    };

    template <typename scalar>
    scalar smooth_heaviside(const scalar &x)
    {
        if (x <= -1)
            return scalar(0.);
        if (x >= 0)
            return scalar(1.);
        return (1. - 2 * x) * intpow(x + 1., 2);
    }

    // template <typename scalar>
    // scalar smooth_heaviside(const scalar &x)
    // {
    //     if (x <= -1)
    //         return scalar(0.);
    //     if (x >= 1)
    //         return scalar(1.);
    //     return (0.5 - x / 4.) * intpow(x + 1., 2);
    // }

    // template <typename scalar>
    // scalar smooth_heaviside(const scalar &x)
    // {
    //     if (x <= 0)
    //         return scalar(0.);
    //     if (x >= 1)
    //         return scalar(1.);
    //     return (3. - 2. * x) * intpow(x, 2);
    // }

    // support is [0, 1]
    template <typename scalar>
    scalar inv_barrier(const scalar &x, const double &r)
    {
        return cubic_spline(2 * x) / pow(x, r);
    }

    template <typename scalar>
    scalar L_ns(const scalar &x)
    {
        if (x < 0.)
            return scalar(0.);
        if (x > 1.)
            return scalar(1.);
        return x;
    }

    template <typename scalar>
    scalar L_s(const scalar &x, const double &a)
    {
        if (x < 0 && x > -a)
            return x - intpow(x / a, 2) * (2 * a + x);
        else
        {
            scalar z = x - 1.;
            if (z < a && z > 0)
                return x - intpow(z / a, 2) * (2 * a - z);
            else
                return L_ns(x);
        }
    }

    template <typename scalar>
    scalar cross2_sqr(const Eigen::Ref<const VectorMax3<scalar>> &a, const Eigen::Ref<const VectorMax3<scalar>> &b)
    {
        assert(a.size() == 3 || a.size() == 2);
        assert(b.size() == 3 || b.size() == 2);
        if (a.size() == 2)
            return intpow(a[0] * b[1] - a[1] * b[0], 2);
        else if (a.size() == 3)
        {
            const Eigen::Ref<const Vector3<scalar>> a_(a), b_(b);
            return a_.cross(b_).squaredNorm();
        }
        else
        {
            assert(false);
            return scalar(0);
        }
    }

    template <typename scalar>
    scalar cross2(const Eigen::Ref<const Vector2<scalar>> &a, const Eigen::Ref<const Vector2<scalar>> &b)
    {
        return a[0] * b[1] - a[1] * b[0];
    }

    // linear solve for 2x2 matrix
    template <typename scalar>
    Vector2<scalar> linear_solve(const Eigen::Ref<const Matrix2<scalar>> &A, const Eigen::Ref<const Vector2<scalar>> &b)
    {
        const scalar det = A(0, 0) * A(1, 1) - A(0, 1) * A(0, 1);
        Vector2<scalar> x;
        x(0) = A(1, 1) * b(0) - A(0, 1) * b(1);
        x(1) = A(1, 0) * b(0) - A(0, 0) * b(1);
        x /= det;
        return x;
    }

    template <class T, int nvert, int dim>
    std::array<Vector<T, dim>, nvert> slice_positions(const Vector<double, nvert*dim> &positions)
    {
        std::array<Vector<T, dim>, nvert> points;
        points.fill(Vector<T, dim>::Zero(dim));
        
        for (int i = 0, id = 0; i < nvert; i++)
            for (int d = 0; d < dim; d++, id++)
                if constexpr (std::is_same<T, double>::value)
                    points[i](d) = positions(id);
                else
                    points[i](d) = T(id, positions(id));

        return points;
    }
}