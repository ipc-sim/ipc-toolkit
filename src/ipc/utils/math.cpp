#include "math.hpp"
#include <ipc/smooth_contact/common.hpp>
#include "AutodiffTypes.hpp"

namespace ipc {
    namespace {
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
                return Math<scalar>::sqr(x) / 2.;
            if (x <= 2)
                return (-3 + x * (6 - 2 * x)) / 2.;
            if (x < 3)
                return Math<scalar>::sqr(3. - x) / 2.;
            return scalar(0.);
        }

        template <typename scalar>
        scalar smooth_heaviside_standard(const scalar &x)
        {
            if (x <= -3)
                return scalar(0.);
            if (x <= -2)
                return Math<scalar>::cubic(3. + x) / 6.;
            if (x <= -1)
                return (((-2 * x - 9) * x - 9) * x + 3) / 6;
            if (x < 0)
                return Math<scalar>::cubic(x) / 6. + 1.;

            return scalar(1.);
        }
    }
    template <typename scalar>
    double Math<scalar>::sign(const double &x)
    {
        if (x > 0)
            return 1.;
        else
            return -1.;
    }

    template <typename scalar>
    scalar Math<scalar>::abs(const scalar &x)
    {
        if (x >= 0)
            return x;
        else
            return -x;
    }

    template <typename scalar>
    scalar Math<scalar>::sqr(const scalar &x)
    {
        return x * x;
    }

    template <typename scalar>
    scalar Math<scalar>::cubic(const scalar &x)
    {
        return x * x * x;
    }

    template <typename scalar>
    scalar Math<scalar>::cubic_spline(const scalar &x)
    {
        if (abs(x) >= 1)
            return scalar(0.);
        if (abs(x) >= 0.5)
            return cubic(1 - abs(x)) * (4. / 3.);
        
        return 2. / 3. - 4. * (x * x) * (1 - abs(x));
    }
    template <typename scalar>
    double Math<scalar>::cubic_spline_grad(const double &x)
    {
        if (Math<double>::abs(x) >= 1)
            return 0.;
        if (Math<double>::abs(x) >= 0.5)
            return -4. * Math<double>::sqr(1 - Math<double>::abs(x)) * sign(x);
        
        return 4. * x * (3. * Math<double>::abs(x) - 2.);
    }
    template <typename scalar>
    double Math<scalar>::cubic_spline_hess(const double &x)
    {
        if (Math<double>::abs(x) >= 1)
            return 0.;
        if (Math<double>::abs(x) >= 0.5)
            return 8. * (1 - Math<double>::abs(x));
        
        return 8. * (3. * Math<double>::abs(x) - 1.);
    }

    /// @brief support is [-1, 1]
    /// @tparam scalar 
    /// @param x 
    /// @return 
    template <typename scalar>
    scalar Math<scalar>::quadratic_spline(const scalar &x)
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
    scalar Math<scalar>::smooth_heaviside(const scalar &x, const double alpha, const double beta)
    {
        return smooth_heaviside_standard((x - beta) * (3 / (alpha + beta)));
    }

    template <typename scalar>
    scalar Math<scalar>::mollifier(const scalar &x)
    {
        if (x <= 0)
            return scalar(0.);
        else if (x < 1)
            return x * (2. - x);
        else
            return scalar(1.);
        // return smooth_heaviside<scalar>(x - 1.);
    }

    template <typename scalar>
    double Math<scalar>::mollifier_grad(const double &x)
    {
        if (x <= 0 || x >= 1)
            return 0.;
        else
            return 2. * (1. - x);
    }

    template <typename scalar>
    double Math<scalar>::mollifier_hess(const double &x)
    {
        if (x <= 0 || x >= 1)
            return 0.;
        else
            return -2.;
    }

    // support is [0, 1]
    template <typename scalar>
    scalar Math<scalar>::inv_barrier(const scalar &x, const double &r)
    {
        return cubic_spline(x) / pow(x, r);
        // log barrier
        // if (x < 1)
        //     return -intpow(1 - sqrt(x), 2) * log(x);
        // else
        //     return scalar(0.);
    }

    template <typename scalar>
    double Math<scalar>::inv_barrier_grad(const double &x, const double &r)
    {
        return (cubic_spline_grad(x) - Math<double>::cubic_spline(x) * r / x) / pow(x, r);
    }
    template <typename scalar>
    double Math<scalar>::inv_barrier_hess(const double &x, const double &r)
    {
        return (cubic_spline_hess(x) + ( -2. * cubic_spline_grad(x) + (r + 1.) * Math<double>::cubic_spline(x) / x ) * r / x) / pow(x, r);
    }

    template <typename scalar>
    scalar Math<scalar>::L_ns(const scalar &x)
    {
        if (x <= 0.)
            return scalar(0.);
        if (x >= 1.)
            return scalar(1.);
        return x;
    }

    template <typename scalar>
    scalar Math<scalar>::cross2(const Eigen::Ref<const Vector2<scalar>> &a, const Eigen::Ref<const Vector2<scalar>> &b)
    {
        return a[0] * b[1] - a[1] * b[0];
    }

    void my_finite_gradient(const Eigen::VectorXd& x, const std::function<double(const Eigen::VectorXd&)> &f, Eigen::VectorXd &grad, FD_RULE rule, const double eps)
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

    template class Math<double>;
    template class Math<ADGrad<1>>;
    template class Math<ADGrad<3>>;
    template class Math<ADGrad<4>>;
    template class Math<ADGrad<6>>;
    template class Math<ADGrad<9>>;
    template class Math<ADGrad<10>>;
    template class Math<ADGrad<12>>;
    template class Math<ADGrad<15>>;
    // template class Math<ADGrad<24>>;
    template class Math<ADGrad<-1, n_vert_neighbors_3d*3+3>>;
    // template class Math<ADGrad<-1, max_vert_3d*3>>;
    template class Math<ADHessian<1>>;
    template class Math<ADHessian<3>>;
    template class Math<ADHessian<4>>;
    template class Math<ADHessian<6>>;
    template class Math<ADHessian<9>>;
    template class Math<ADHessian<10>>;
    template class Math<ADHessian<12>>;
    template class Math<ADHessian<15>>;
    // template class Math<ADHessian<24>>;
    template class Math<ADHessian<-1, n_vert_neighbors_3d*3+3>>;
    // template class Math<ADHessian<-1, max_vert_3d*3>>;
}