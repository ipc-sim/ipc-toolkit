#include "edge_edge.hpp"

#include <unsupported/Eigen/KroneckerProduct>
#include <ipc/friction/closest_point.hpp>
#include "autogen.hpp"

namespace ipc {

std::tuple<Vector3d, Eigen::Matrix<double, 3, 12>>
line_line_closest_point_direction_gradient(
    const Eigen::Ref<const Vector3d>& ea0,
    const Eigen::Ref<const Vector3d>& ea1,
    const Eigen::Ref<const Vector3d>& eb0,
    const Eigen::Ref<const Vector3d>& eb1)
{
    const Vector3d ea = ea1 - ea0;
    const Vector3d eb = eb1 - eb0;
    const Vector3d normal = ea.cross(eb);
    const Eigen::Matrix<double, 3, 6> cross_grad = cross_product_gradient(ea, eb);

    const double normal_sqr_norm = normal.squaredNorm();
    const Vector3d t = eb0 - ea0;
    const Vector3d vec = (t.dot(normal) / normal_sqr_norm) * normal;

    Eigen::Matrix<double, 3, 12> grad = Eigen::Matrix<double, 3, 12>::Zero();
    // derivative wrt. normal, tempararily store here
    grad(Eigen::all, {3,4,5}) = (normal * t.transpose() + t.transpose() * normal * Eigen::Matrix3d::Identity() - 2 * vec * normal.transpose()) / normal_sqr_norm;
    // derivative wrt. t
    grad.middleCols<3>(6) = normal * normal.transpose() / normal_sqr_norm;
    grad.leftCols<3>() = -grad.middleCols<3>(6);
    // contributions from n
    grad(Eigen::all, {3,4,5, 9,10,11}) = grad(Eigen::all, {3,4,5}) * cross_grad;
    grad(Eigen::all, {0,1,2, 6,7,8})   -= grad(Eigen::all, {3,4,5, 9,10,11});

    return std::make_tuple(vec, grad);
}

std::tuple<Vector3d, Eigen::Matrix<double, 3, 12>, std::array<Matrix12d, 3>>
line_line_closest_point_direction_hessian(
    const Eigen::Ref<const Vector3d>& ea0,
    const Eigen::Ref<const Vector3d>& ea1,
    const Eigen::Ref<const Vector3d>& eb0,
    const Eigen::Ref<const Vector3d>& eb1)
{
    const Vector3d ea = ea1 - ea0;
    const Vector3d eb = eb1 - eb0;
    const Vector3d normal = ea.cross(eb);
    const Eigen::Matrix<double, 3, 6> cross_grad = cross_product_gradient(ea, eb);
    const std::array<Matrix6d, 3> cross_hess = cross_product_hessian(ea, eb);

    const double normal_sqr_norm = normal.squaredNorm();
    const Vector3d t = eb0 - ea0;
    const Vector3d vec = (t.dot(normal) / normal_sqr_norm) * normal;

    Eigen::Matrix<double, 3, 12> grad = Eigen::Matrix<double, 3, 12>::Zero();
    // derivative wrt. normal and t
    Eigen::Matrix3d grad_normal = (normal * t.transpose() - 2 * vec * normal.transpose()) / normal_sqr_norm;
    grad_normal.diagonal().array() += (t.transpose() * normal)(0) / normal_sqr_norm;
    Eigen::Matrix3d grad_t = normal * normal.transpose() / normal_sqr_norm;
    grad.middleCols<3>(6) = grad_t;
    grad.leftCols<3>() = -grad_t;
    // contributions from n to ea, eb
    grad(Eigen::all, {3,4,5, 9,10,11}) = grad_normal * cross_grad;
    grad(Eigen::all, {0,1,2, 6,7,8})   -= grad(Eigen::all, {3,4,5, 9,10,11});

    // hessian of output wrt. normal
    std::array<Eigen::Matrix3d, 3> hess_wrt_normal;
    // hessian of output wrt. normal (row) and t (col)
    std::array<Eigen::Matrix3d, 3> mixed_hess;
    for (int d = 0; d < 3; d++)
    {
        hess_wrt_normal[d] = -2 * (normal * grad_normal.row(d) + grad_normal.row(d).transpose() * normal.transpose());
        hess_wrt_normal[d].col(d) += t;
        hess_wrt_normal[d].row(d) += t.transpose();
        hess_wrt_normal[d].diagonal().array() -= 2 * vec(d);
        hess_wrt_normal[d] /= normal_sqr_norm;

        mixed_hess[d] = -2 * normal * grad_t.row(d);
        mixed_hess[d].row(d) += normal.transpose();
        mixed_hess[d].diagonal().array() += normal(d);
        mixed_hess[d] /= normal_sqr_norm;
    }

    // grad and hessian of [normal, t] wrt. [ea0, ea1, eb0, eb1]
    Eigen::Matrix<double, 3, 12> inner_grad;
    inner_grad << -cross_grad.leftCols<3>(), cross_grad.leftCols<3>(), -cross_grad.rightCols<3>(), cross_grad.rightCols<3>();
    
    std::array<Matrix12d, 3> inner_hess;
    {
        Eigen::Matrix2d tmp;
        tmp << 1, -1, -1, 1;
        for (int d = 0; d < 3; d++)
            inner_hess[d] << Eigen::KroneckerProduct<Eigen::Matrix2d, Eigen::Matrix3d>(tmp, cross_hess[d].block<3, 3>(0, 0)),
                        Eigen::KroneckerProduct<Eigen::Matrix2d, Eigen::Matrix3d>(tmp, cross_hess[d].block<3, 3>(0, 3)),
                        Eigen::KroneckerProduct<Eigen::Matrix2d, Eigen::Matrix3d>(tmp, cross_hess[d].block<3, 3>(3, 0)),
                        Eigen::KroneckerProduct<Eigen::Matrix2d, Eigen::Matrix3d>(tmp, cross_hess[d].block<3, 3>(3, 3));
    }
    
    std::array<Matrix12d, 3> hess;
    for (int d = 0; d < 3; d++)
    {
        Eigen::Matrix<double, 12, 3> tmp = inner_grad.transpose() * mixed_hess[d];
        hess[d] = inner_grad.transpose() * hess_wrt_normal[d] * inner_grad;
        hess[d].middleCols<3>(6) += tmp;
        hess[d].middleCols<3>(0) -= tmp;
        hess[d].middleRows<3>(6) += tmp.transpose();
        hess[d].middleRows<3>(0) -= tmp.transpose();
        for (int i = 0; i < 3; i++)
            hess[d] += inner_hess[i] * grad_normal(d, i);
    }

    return std::make_tuple(vec, grad, hess);
}


std::tuple<Vector6d, Eigen::Matrix<double, 6, 12>>
line_line_closest_point_pairs_gradient(
    const Eigen::Ref<const Vector3d>& ea0,
    const Eigen::Ref<const Vector3d>& ea1,
    const Eigen::Ref<const Vector3d>& eb0,
    const Eigen::Ref<const Vector3d>& eb1)
{
    const auto uv = edge_edge_closest_point(ea0, ea1, eb0, eb1);
    const auto J = edge_edge_closest_point_jacobian(ea0, ea1, eb0, eb1);

    Vector6d out;
    out << uv(0) * (ea1 - ea0) + ea0, uv(1) * (eb1 - eb0) + eb0;
    Eigen::Matrix<double, 6, 12> grad = Eigen::Matrix<double, 6, 12>::Zero();
    grad.topRows<3>() += (ea1 - ea0) * J.row(0);
    grad.block<3, 3>(0, 3).diagonal().array() += uv(0);
    grad.block<3, 3>(0, 0).diagonal().array() += (1 - uv(0));
    grad.bottomRows<3>() += (eb1 - eb0) * J.row(1);
    grad.block<3, 3>(3, 9).diagonal().array() += uv(1);
    grad.block<3, 3>(3, 6).diagonal().array() += (1 - uv(1));

    return {out, grad};
}

std::tuple<Vector6d, Eigen::Matrix<double, 6, 12>, std::array<Matrix12d, 6>>
line_line_closest_point_pairs_hessian(
    const Eigen::Ref<const Vector3d>& ea0,
    const Eigen::Ref<const Vector3d>& ea1,
    const Eigen::Ref<const Vector3d>& eb0,
    const Eigen::Ref<const Vector3d>& eb1)
{
    const auto uv = edge_edge_closest_point(ea0, ea1, eb0, eb1);
    const auto J = edge_edge_closest_point_jacobian(ea0, ea1, eb0, eb1);

    Vector6d out;
    out << uv(0) * (ea1 - ea0) + ea0, uv(1) * (eb1 - eb0) + eb0;

    Eigen::Matrix<double, 6, 12> grad = Eigen::Matrix<double, 6, 12>::Zero();
    grad.topRows<3>() += (ea1 - ea0) * J.row(0);
    grad.block<3, 3>(0, 3).diagonal().array() += uv(0);
    grad.block<3, 3>(0, 0).diagonal().array() += (1 - uv(0));
    grad.bottomRows<3>() += (eb1 - eb0) * J.row(1);
    grad.block<3, 3>(3, 9).diagonal().array() += uv(1);
    grad.block<3, 3>(3, 6).diagonal().array() += (1 - uv(1));

    std::array<Matrix12d, 6> hess;
    for (auto &h : hess)
        h.setZero();
    {
        Eigen::Matrix<double, 12, 12> Ha, Hb;
        autogen::edge_edge_closest_point_hessian_a(
            ea0.x(), ea0.y(), ea0.z(),
            ea1.x(), ea1.y(), ea1.z(),
            eb0.x(), eb0.y(), eb0.z(),
            eb1.x(), eb1.y(), eb1.z(),
            Ha.data());
        autogen::edge_edge_closest_point_hessian_b(
            ea0.x(), ea0.y(), ea0.z(),
            ea1.x(), ea1.y(), ea1.z(),
            eb0.x(), eb0.y(), eb0.z(),
            eb1.x(), eb1.y(), eb1.z(),
            Hb.data());
        
        for (int d = 0; d < 3; d++)
        {
            // wrt. uv
            hess[d] += (ea1(d) - ea0(d)) * Ha;
            hess[d+3] += (eb1(d) - eb0(d)) * Hb;

            // wrt. ea0 & uv(0)
            hess[d].row(d) -= J.row(0);
            hess[d].col(d) -= J.row(0).transpose();

            // wrt. ea1 & uv(0)
            hess[d].row(d+3) += J.row(0);
            hess[d].col(d+3) += J.row(0).transpose();

            // wrt. eb0 & uv(1)
            hess[d+3].row(d+6) -= J.row(1);
            hess[d+3].col(d+6) -= J.row(1).transpose();

            // wrt. eb1 & uv(1)
            hess[d+3].row(d+9) += J.row(1);
            hess[d+3].col(d+9) += J.row(1).transpose();
        }
    }

    return {out, grad, hess};
}
}
