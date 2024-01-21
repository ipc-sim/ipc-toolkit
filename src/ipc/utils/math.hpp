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
        scalar y = 2. * x;
        if (y <= -2)
            return scalar(0.);
        if (y <= -1)
            return intpow(y + 2, 3) / 6;
        if (y <= 0)
            return 2. / 3 - intpow(y, 2) * (1 + y / 2);
        if (y <= 1)
            return 2. / 3 - intpow(y, 2) * (1 - y / 2);
        if (y < 2)
            return -intpow(y - 2, 3) / 6;
        return scalar(0.);
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
            return intpow(x, 2) / 2.;
        if (x <= 2)
            return (-3 + x * (6 - 2 * x)) / 2.;
        if (x < 3)
            return intpow(3. - x, 2) / 2.;
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

    enum class HEAVISIDE_TYPE {
        ZERO = 0, ONE = 1, VARIANT = 2
    };

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
    scalar smooth_heaviside_aux(const scalar &x)
    {
        if (x <= -3)
            return scalar(0.);
        if (x <= -2)
            return intpow(3. + x, 3) / 6.;
        if (x <= -1)
            return (((-2 * x - 9) * x - 9) * x + 3) / 6;
        if (x < 0)
            return intpow(x, 3) / 6. + 1.;

        return scalar(1.);
    }

    template <typename scalar>
    scalar smooth_heaviside(const scalar &x)
    {
        return smooth_heaviside_aux(3 * x);
    }

    constexpr double mollifier_threshold_eps = 1e-3;

    template <typename scalar>
    scalar mollifier(const scalar &x)
    {
        // if (x <= 0)
        //     return scalar(0.);
        // if (x <= 1)
        //     return x * (2. - x);
        // return scalar(1.);
        return smooth_heaviside<scalar>(x - 1.);
    }

    // support is [0, 1]
    template <typename scalar>
    scalar inv_barrier(const scalar &x, const double &r)
    {
        return cubic_spline<scalar>(x) / pow(x, r);
        // if (x < 1)
        //     return -pow(1 - x, 2.) * log(x);
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

    template <typename scalar>
    scalar normalized_cross2(const Eigen::Ref<const Vector2<scalar>> &a, const Eigen::Ref<const Vector2<scalar>> &b)
    {
        return cross2(a, b) / a.norm() / b.norm();
    }

    template <typename scalar, int dim>
    scalar normalized_dot(const Eigen::Ref<const Vector<scalar, dim>> &a, const Eigen::Ref<const Vector<scalar, dim>> &b)
    {
        return a.dot(b) / a.norm() / b.norm();
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

    template <class T, int dim>
    Eigen::Matrix<T, -1, dim> slice_positions_large(const Eigen::VectorXd &positions)
    {
        const int nvert = positions.size() / dim;
        Eigen::Matrix<T, -1, dim> points;
        points.setZero(nvert, dim);
        
        for (int i = 0, id = 0; i < nvert; i++)
            for (int d = 0; d < dim; d++, id++)
                if constexpr (std::is_same<T, double>::value)
                    points(i, d) = positions(id);
                else
                    points(i, d) = T(id, positions(id));

        return points;
    }

    template <typename Less, typename T, typename... Ts>
    inline constexpr const T& min(Less less, const T& a, const T& b, const Ts&... rems) {
        if constexpr (sizeof...(rems)) {
            return min(less, std::min(a, b, less), rems...);
        }
        else {
            return std::min(a, b, less);
        }
    }

    template <typename Less, typename T, typename... Ts>
    inline constexpr const T& max(Less less, const T& a, const T& b, const Ts&... rems) {
        if constexpr (sizeof...(rems)) {
            return max(less, std::max(a, b, less), rems...);
        }
        else {
            return std::max(a, b, less);
        }
    }

    enum class FD_RULE { CENTRAL, LEFT, RIGHT };
    
    inline void my_finite_gradient(const Eigen::VectorXd& x, const std::function<double(const Eigen::VectorXd&)> &f, Eigen::VectorXd &grad, FD_RULE rule = FD_RULE::CENTRAL, const double eps = 1e-7)
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