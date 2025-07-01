#include "math.hpp"

#include "AutodiffTypes.hpp"

#include <ipc/smooth_contact/common.hpp>

DECLARE_DIFFSCALAR_BASE();

namespace ipc {
namespace {

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

GradType<6> opposite_direction_penalty_grad(
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

HessianType<6> opposite_direction_penalty_hess(
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

GradType<9> negative_orientation_penalty_grad(
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

HessianType<9> negative_orientation_penalty_hess(
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2,
    const Eigen::Ref<const Eigen::Vector3d>& d,
    const double& alpha,
    const double& beta)
{
    const Vector3d n = t1.cross(t2);
    const Eigen::Matrix<double, 3, 6> cross_grad =
        cross_product_gradient(t1, t2);
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

} // namespace ipc