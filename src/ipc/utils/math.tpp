#pragma once
#include "math.hpp"

namespace ipc {
    template <typename scalar>
    scalar abs(const scalar &x)
    {
        if (x >= 0)
            return x;
        else
            return -x;
    }

    template <typename scalar>
    scalar cubic_spline(const scalar &x)
    {
        const scalar y = 2. * x;
        if (abs(y) >= 2)
            return scalar(0.);
        if (abs(y) >= 1)
            return cubic(2 - abs(y)) / 6;
        
        return 2. / 3 - (y * y) * (1 - abs(y) / 2);
    }

    /// @brief support is [0, 3]
    /// @tparam scalar 
    /// @param x 
    /// @return 
    template <typename scalar>
    scalar quadratic_spline_aux(const scalar &x)
    {
        if (x <= 0)
            return scalar(0.);
        if (x <= 1)
            return x * x / 2.;
        if (x <= 2)
            return (-3 + x * (6 - 2 * x)) / 2.;
        if (x < 3)
            return sqr(3. - x) / 2.;
        return scalar(0.);
    }

    /// @brief support is [-1, 1]
    /// @tparam scalar 
    /// @param x 
    /// @return 
    template <typename scalar>
    scalar quadratic_spline(const scalar &x)
    {
        return quadratic_spline_aux(x * 1.5 + 1.5);
    }

    HEAVISIDE_TYPE ORIENTATION_TYPES::compute_type(const double &val, const double &alpha, const double &beta)
    {
        if (val <= -alpha)
            return HEAVISIDE_TYPE::ZERO;
        if (val >= beta)
            return HEAVISIDE_TYPE::ONE;
        
        return HEAVISIDE_TYPE::VARIANT;
    }

    void ORIENTATION_TYPES::set_size(const int size)
    {
        size_ = size;
        tangent_types.assign(size_, HEAVISIDE_TYPE::VARIANT);
        normal_types.assign(size_, HEAVISIDE_TYPE::VARIANT);
    }

    // template <typename scalar>
    // scalar smooth_heaviside_aux(const scalar &x)
    // {
    //     if (x <= -2)
    //         return scalar(0.);
    //     if (x <= -1)
    //         return intpow(2. + x, 2) / 2.;
    //     if (x <= 0)
    //         return 1 - intpow(x, 2) / 2.;

    //     return scalar(1.);
    // }

    template <typename scalar>
    scalar smooth_heaviside_standard(const scalar &x)
    {
        if (x <= -3)
            return scalar(0.);
        if (x <= -2)
            return cubic(3. + x) / 6.;
        if (x <= -1)
            return (((-2 * x - 9) * x - 9) * x + 3) / 6;
        if (x < 0)
            return cubic(x) / 6. + 1.;

        return scalar(1.);
    }

    template <typename scalar>
    scalar smooth_heaviside(const scalar &x, const double alpha, const double beta)
    {
        return smooth_heaviside_standard((x - beta) * (3 / (alpha + beta)));
    }

    template <typename scalar>
    scalar mollifier(const scalar &x)
    {
        if (x <= 0)
            return scalar(0.);
        if (x <= 1)
        {
            return x * (2. - x);
        }
        return scalar(1.);
        // return smooth_heaviside<scalar>(x - 1., 1., 0.);
    }

    // support is [0, 1]
    template <typename scalar>
    scalar inv_barrier(const scalar &x, const double &r)
    {
        return cubic_spline<scalar>(x) / pow(x, r);
        // if (x < 1)
        //     return -intpow(1 - sqrt(x), 2) * log(x);
        // else
        //     return scalar(0.);
    }

    template <typename scalar>
    scalar L_ns(const scalar &x)
    {
        if (x <= 0.)
            return scalar(0.);
        if (x >= 1.)
            return scalar(1.);
        return x;
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

    template <class T, int rows, int dim, int max_rows=rows>
    Eigen::Matrix<T, rows, dim, Eigen::ColMajor, max_rows, dim> slice_positions(const Eigen::VectorXd &positions)
    {
        const int nvert = positions.size() / dim;
        Eigen::Matrix<T, rows, dim, Eigen::ColMajor, max_rows, dim> points;
        points.setZero(nvert, dim);
        
        for (int i = 0, id = 0; i < nvert; i++)
            for (int d = 0; d < dim; d++, id++)
                if constexpr (std::is_same<T, double>::value)
                    points(i, d) = positions(id);
                else
                    points(i, d) = T(id, positions(id));

        return points;
    }

    inline void my_finite_gradient(const Eigen::VectorXd& x, const std::function<double(const Eigen::VectorXd&)> &f, Eigen::VectorXd &grad, FD_RULE rule, const double eps)
    {
        grad.setZero(x.size());
        switch (rule)
        {
        case FD_RULE::CENTRAL:
            for (int i = 0; i < x.size(); i++)
                for (int d : {-1, 1})
                {
                    auto y = x;
                    y(i) += d * eps;
                    grad(i) += d * f(y) / (2*eps);
                }
            break;
        case FD_RULE::LEFT:
            for (int i = 0; i < x.size(); i++)
            {
                    auto y = x;
                    grad(i) += f(y) / eps;
                    y(i) -= eps;
                    grad(i) -= f(y) / eps;
            }
            break;
        case FD_RULE::RIGHT:
            for (int i = 0; i < x.size(); i++)
            {
                    auto y = x;
                    grad(i) -= f(y) / eps;
                    y(i) += eps;
                    grad(i) += f(y) / eps;
            }
            break;
        default:
        assert(false);
        }
    }
}