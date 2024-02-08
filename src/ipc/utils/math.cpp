#include "math.hpp"
#include <ipc/smooth_contact/common.hpp>
#include "AutodiffTypes.hpp"

namespace ipc {
namespace {
    /// @brief support is [0, 3]
    /// @tparam scalar
    /// @param x
    /// @return
    template <typename scalar> scalar quadratic_spline_aux(const scalar& x)
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

    template <typename scalar> scalar smooth_heaviside_standard(const scalar& x)
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

    double smooth_heaviside_standard_grad(const double& x)
    {
        if (x <= -3 || x >= 0)
            return 0.;
        if (x <= -2)
            return Math<double>::sqr(3. + x) / 2.;
        if (x <= -1)
            return -(x * x + 3 * x + 1.5);

        return Math<double>::sqr(x) / 2.;
    }

    double smooth_heaviside_standard_hess(const double& x)
    {
        if (x <= -3 || x >= 0)
            return 0.;
        if (x <= -2)
            return 3. + x;
        if (x <= -1)
            return -3 - 2 * x;

        return x;
    }

    /*
        Mathematica script

        t1x = Table[t1[k], {k, 0, 2}];
        t2x = Table[t2[k], {k, 0, 2}];
        t = Flatten[{t1x, t2x}];
        n = Cross[t1x, t2x];
        CForm[Sum[ddy[i - 1, j - 1] Transpose[{D[n[[j]], {t}]}].{D[n[[i]], {t}]}
       + dy[i - 1] D[n[[i]], {t}, {t}], {i, 1, 3}, {j, 1, 3}]]
    */

} // namespace
template <typename scalar> double Math<scalar>::sign(const double& x)
{
    if (x > 0)
        return 1.;
    else
        return -1.;
}

template <typename scalar> scalar Math<scalar>::abs(const scalar& x)
{
    if (x >= 0)
        return x;
    else
        return -x;
}

template <typename scalar> scalar Math<scalar>::sqr(const scalar& x)
{
    return x * x;
}

template <typename scalar> scalar Math<scalar>::cubic(const scalar& x)
{
    return x * x * x;
}

template <typename scalar> scalar Math<scalar>::cubic_spline(const scalar& x)
{
    if (abs(x) >= 1)
        return scalar(0.);
    if (abs(x) >= 0.5)
        return cubic(1 - abs(x)) * (4. / 3.);

    return 2. / 3. - 4. * (x * x) * (1 - abs(x));
}
template <typename scalar>
double Math<scalar>::cubic_spline_grad(const double& x)
{
    if (Math<double>::abs(x) >= 1)
        return 0.;
    if (Math<double>::abs(x) >= 0.5)
        return -4. * Math<double>::sqr(1 - Math<double>::abs(x)) * sign(x);

    return 4. * x * (3. * Math<double>::abs(x) - 2.);
}
template <typename scalar>
double Math<scalar>::cubic_spline_hess(const double& x)
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
scalar Math<scalar>::quadratic_spline(const scalar& x)
{
    return quadratic_spline_aux(x * 1.5 + 1.5);
}

HEAVISIDE_TYPE ORIENTATION_TYPES::compute_type(
    const double& val, const double& alpha, const double& beta)
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

bool ORIENTATION_TYPES::are_tangent_types_all_one() const
{
    for (const auto& b : tangent_types)
        if (b != HEAVISIDE_TYPE::ONE)
            return false;
    return true;
}

bool ORIENTATION_TYPES::exists_normal_type_one() const
{
    for (const auto& b : normal_types)
        if (b == HEAVISIDE_TYPE::ONE)
            return true;
    return false;
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
scalar Math<scalar>::smooth_heaviside(
    const scalar& x, const double alpha, const double beta)
{
    return smooth_heaviside_standard((x - beta) * (3 / (alpha + beta)));
}

template <typename scalar>
double Math<scalar>::smooth_heaviside_grad(
    const double& x, const double alpha, const double beta)
{
    const double s = 3 / (alpha + beta);
    return smooth_heaviside_standard_grad((x - beta) * s) * s;
}

template <typename scalar>
double Math<scalar>::smooth_heaviside_hess(
    const double& x, const double alpha, const double beta)
{
    const double s = 3 / (alpha + beta);
    return smooth_heaviside_standard_hess((x - beta) * s) * s * s;
}

template <typename scalar> scalar Math<scalar>::mollifier(const scalar& x)
{
    if constexpr (isADHessian<scalar>::value) {
        if (x <= 0)
            return scalar(0.);
        else if (x < 1) {
            const double deriv = 2. * (1. - x.getValue()), hess = -2.;
            return scalar(
                x.getValue() * (2. - x.getValue()), deriv * x.getGradient(),
                x.getGradient() * hess * x.getGradient().transpose()
                    + deriv * x.getHessian());
        } else
            return scalar(1.);
        // return smooth_heaviside<scalar>(x - 1.);
    } else {
        if (x <= 0)
            return scalar(0.);
        else if (x < 1)
            return x * (2. - x);
        else
            return scalar(1.);
        // return smooth_heaviside<scalar>(x - 1.);
    }
}

template <typename scalar> double Math<scalar>::mollifier_grad(const double& x)
{
    if (x <= 0 || x >= 1)
        return 0.;
    else
        return 2. * (1. - x);
}

template <typename scalar> double Math<scalar>::mollifier_hess(const double& x)
{
    if (x <= 0 || x >= 1)
        return 0.;
    else
        return -2.;
}

// support is [0, 1]
template <typename scalar>
scalar Math<scalar>::inv_barrier(const scalar& x, const double& r)
{
    return cubic_spline(x) / pow(x, r);
    // log barrier
    // if (x < 1)
    //     return -intpow(1 - sqrt(x), 2) * log(x);
    // else
    //     return scalar(0.);
}

template <typename scalar>
double Math<scalar>::inv_barrier_grad(const double& x, const double& r)
{
    return (cubic_spline_grad(x) - Math<double>::cubic_spline(x) * r / x)
        / pow(x, r);
}
template <typename scalar>
double Math<scalar>::inv_barrier_hess(const double& x, const double& r)
{
    return (cubic_spline_hess(x)
            + (-2. * cubic_spline_grad(x)
               + (r + 1.) * Math<double>::cubic_spline(x) / x)
                * r / x)
        / pow(x, r);
}

template <typename scalar> scalar Math<scalar>::L_ns(const scalar& x)
{
    if (x <= 0.)
        return scalar(0.);
    if (x >= 1.)
        return scalar(1.);
    return x;
}

template <typename scalar>
scalar Math<scalar>::cross2(
    const Eigen::Ref<const Vector2<scalar>>& a,
    const Eigen::Ref<const Vector2<scalar>>& b)
{
    return a[0] * b[1] - a[1] * b[0];
}

void my_finite_gradient(
    const Eigen::VectorXd& x,
    const std::function<double(const Eigen::VectorXd&)>& f,
    Eigen::VectorXd& grad,
    FD_RULE rule,
    const double eps)
{
    grad.setZero(x.size());
    switch (rule) {
    case FD_RULE::CENTRAL:
        for (int i = 0; i < x.size(); i++)
            for (int d : { -1, 1 }) {
                auto y = x;
                y(i) += d * eps;
                grad(i) += d * f(y) / (2 * eps);
            }
        break;
    case FD_RULE::LEFT:
        for (int i = 0; i < x.size(); i++) {
            auto y = x;
            grad(i) += f(y) / eps;
            y(i) -= eps;
            grad(i) -= f(y) / eps;
        }
        break;
    case FD_RULE::RIGHT:
        for (int i = 0; i < x.size(); i++) {
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

std::tuple<Eigen::Vector3d, Eigen::Matrix3d>
normalize_vector_grad(const Eigen::Ref<const Eigen::Vector3d>& t)
{
    double norm = t.norm();
    Eigen::Vector3d y = t / norm;
    Eigen::Matrix3d grad =
        (Eigen::Matrix3d::Identity() - y * y.transpose()) / norm;
    return std::make_tuple(y, grad);
}

std::tuple<
    Eigen::Vector3d,
    Eigen::Matrix3d,
    std::array<Eigen::Matrix<double, 3, 3>, 3>>
normalize_vector_hess(const Eigen::Ref<const Eigen::Vector3d>& t)
{
    double norm = t.norm();
    Eigen::Vector3d y = t / norm;
    Eigen::Matrix3d grad =
        (Eigen::Matrix3d::Identity() - y * y.transpose()) / norm;
    std::array<Eigen::Matrix<double, 3, 3>, 3> hess;
    for (int i = 0; i < 3; i++)
        hess[i] = -(y(i) * grad + y * grad.row(i) + grad.col(i) * y.transpose())
            / norm;

    return std::make_tuple(y, grad, hess);
}

double opposite_direction_penalty(
    const Eigen::Ref<const Eigen::Vector3d>& t,
    const Eigen::Ref<const Eigen::Vector3d>& d,
    const double& alpha,
    const double& beta)
{
    return Math<double>::smooth_heaviside(d.dot(t) / t.norm(), alpha, beta);
}

std::tuple<double, Vector6d> opposite_direction_penalty_grad(
    const Eigen::Ref<const Eigen::Vector3d>& t,
    const Eigen::Ref<const Eigen::Vector3d>& d,
    const double& alpha,
    const double& beta)
{
    auto [tn, tn_grad] = normalize_vector_grad(t);
    const double a = d.dot(tn);
    const double y = Math<double>::smooth_heaviside(a, alpha, beta);
    const double dy = Math<double>::smooth_heaviside_grad(a, alpha, beta);

    Vector6d grad;
    grad << tn_grad * dy * d, dy * tn;
    return std::make_tuple(y, grad);
}

std::tuple<double, Vector6d, Matrix6d> opposite_direction_penalty_hess(
    const Eigen::Ref<const Eigen::Vector3d>& t,
    const Eigen::Ref<const Eigen::Vector3d>& d,
    const double& alpha,
    const double& beta)
{
    auto [tn, tn_grad, tn_hess] = normalize_vector_hess(t);
    const double a = d.dot(tn);
    const double y = Math<double>::smooth_heaviside(a, alpha, beta);
    const double dy = Math<double>::smooth_heaviside_grad(a, alpha, beta);
    const double ddy = Math<double>::smooth_heaviside_hess(a, alpha, beta);

    Vector6d grad;
    grad << dy * d, dy * tn;

    Matrix6d hess;
    hess << d * ddy * d.transpose(),
        d * ddy * tn.transpose() + Eigen::Matrix3d::Identity() * dy,
        tn * ddy * d.transpose() + Eigen::Matrix3d::Identity() * dy,
        tn * ddy * tn.transpose();

    // chain rule of vector t normalize
    hess.topRows(3) = tn_grad * hess.topRows(3);
    hess.leftCols(3) = hess.leftCols(3) * tn_grad;
    hess.topLeftCorner(3, 3) +=
        grad(0) * tn_hess[0] + grad(1) * tn_hess[1] + grad(2) * tn_hess[2];

    grad.head(3) = tn_grad * grad.head(3);

    return std::make_tuple(y, grad, hess);
}

// assume unit vector d
double negative_orientation_penalty(
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2,
    const Eigen::Ref<const Eigen::Vector3d>& d,
    const double& alpha,
    const double& beta)
{
    const Eigen::Vector3d n = t1.cross(t2);
    return opposite_direction_penalty(n, d, alpha, beta);
}

std::tuple<double, Vector9d> negative_orientation_penalty_grad(
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2,
    const Eigen::Ref<const Eigen::Vector3d>& d,
    const double& alpha,
    const double& beta)
{
    const Eigen::Vector3d n = t1.cross(t2);
    auto [y, dy] = opposite_direction_penalty_grad(n, d, alpha, beta);

    /*
        Mathematica script

        t1 = {t10, t11, t12};
        t2 = {t20, t21, t22};
        n = Cross[t1, t2];
        Sum[dy[k - 1] D[n[[k]], {{t1, t2}}], {k, 1, 3}] // MatrixForm
    */

    Vector9d grad;
    grad << -t2(2) * dy(1) + t2(1) * dy(2), -t2(0) * dy(2) + t2(2) * dy(0),
        -t2(1) * dy(0) + t2(0) * dy(1), -t1(1) * dy(2) + t1(2) * dy(1),
        -t1(2) * dy(0) + t1(0) * dy(2), -t1(0) * dy(1) + t1(1) * dy(0),
        dy.tail(3);
    return std::make_tuple(y, grad);
}

std::tuple<double, Vector9d, Matrix9d> negative_orientation_penalty_hess(
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2,
    const Eigen::Ref<const Eigen::Vector3d>& d,
    const double& alpha,
    const double& beta)
{
    const Vector3d n = t1.cross(t2);
    const Eigen::Matrix<double, 3, 6> cross_grad = cross_product_gradient(t1, t2);
    const std::array<Matrix6d, 3> cross_hess = cross_product_hessian(t1, t2);
    auto [y, dy, ddy] = opposite_direction_penalty_hess(n, d, alpha, beta);

    Vector9d grad;
    grad.tail(3) = dy.tail(3);
    grad.head(6) = cross_grad.transpose() * dy.head(3);

    Matrix9d hess;
    hess.bottomRightCorner(3, 3) = ddy.bottomRightCorner(3, 3);
    hess.topLeftCorner(6, 6) =
        cross_grad.transpose() * ddy.topLeftCorner(3, 3) * cross_grad
        + cross_hess[0] * dy(0) + cross_hess[1] * dy(1) + cross_hess[2] * dy(2);
    hess.topRightCorner(6, 3) =
        cross_grad.transpose() * ddy.topRightCorner(3, 3);
    hess.bottomLeftCorner(3, 6) = ddy.bottomLeftCorner(3, 3) * cross_grad;

    return std::make_tuple(y, grad, hess);
}

Eigen::Matrix<double, 3, 6> cross_product_gradient(
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2)
{
    Eigen::Matrix<double, 3, 6> grad;
    grad << 0, t2(2), -t2(1), 0, -t1(2), t1(1), -t2(2), 0, t2(0), t1(2), 0,
        -t1(0), t2(1), -t2(0), 0, -t1(1), t1(0), 0;

    return grad;
}

std::array<Matrix6d, 3> cross_product_hessian(
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2)
{
    std::array<Matrix6d, 3> hess;
    hess.fill(Matrix6d::Zero());
    hess[0](1, 5) = 1;
    hess[0](5, 1) = 1;
    hess[0](2, 4) = -1;
    hess[0](4, 2) = -1;

    hess[1](0, 5) = -1;
    hess[1](5, 0) = -1;
    hess[1](2, 3) = 1;
    hess[1](3, 2) = 1;

    hess[2](0, 4) = 1;
    hess[2](4, 0) = 1;
    hess[2](1, 3) = -1;
    hess[2](3, 1) = -1;

    return hess;
}


template class Math<double>;

template class Math<ADGrad<1>>;
template class Math<ADHessian<1>>;

template class Math<ADGrad<3>>;
template class Math<ADHessian<3>>;

template class Math<ADGrad<4>>;
template class Math<ADHessian<4>>;

template class Math<ADGrad<6>>;
template class Math<ADHessian<6>>;

template class Math<ADGrad<8>>;
template class Math<ADHessian<8>>;

template class Math<ADGrad<9>>;
template class Math<ADHessian<9>>;

template class Math<ADGrad<10>>;
template class Math<ADHessian<10>>;

template class Math<ADGrad<12>>;
template class Math<ADHessian<12>>;

// #ifdef DERIVATIVES_WITH_AUTODIFF
template class Math<ADGrad<13>>;
template class Math<ADHessian<13>>;

template class Math<ADGrad<15>>;
template class Math<ADHessian<15>>;

template class Math<ADGrad<-1>>;
template class Math<ADHessian<-1>>;
// #endif
} // namespace ipc