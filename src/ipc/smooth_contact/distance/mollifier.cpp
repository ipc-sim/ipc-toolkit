#ifndef DERIVATIVES_WITH_AUTODIFF

#include "mollifier.hpp"

#include "point_edge.hpp"

#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_point.hpp>

namespace ipc {
namespace {
    // double func_aux(const double& a, const double& b, const double& c, const
    // double& eps)
    // {
    //     return Math<double>::mollifier((a - c) / c / eps);
    // }

    Vector<int, 9> get_indices(int i, int j, int k)
    {
        Vector<int, 9> out;
        out << 3 * i + 0, 3 * i + 1, 3 * i + 2, 3 * j + 0, 3 * j + 1, 3 * j + 2,
            3 * k + 0, 3 * k + 1, 3 * k + 2;
        return out;
    }

    GradType<3> func_aux_grad(
        const double& a, const double& b, const double& c, const double& eps)
    {
        const double val = (a - c) / c / eps;
        return std::make_tuple(
            Math<double>::mollifier(val),
            Eigen::Vector3d(1. / c, 0., -a / c / c)
                * Math<double>::mollifier_grad(val) / eps);
    }

    HessianType<3> func_aux_hess(
        const double& a, const double& b, const double& c, const double& eps)
    {
        const double c2 = 1. / c / c;
        Eigen::Matrix3d h1;
        h1 << 0., 0., -c2, 0., 0., 0., -c2, 0., 2 * a * c2 / c;
        Eigen::Vector3d g1;
        g1 << 1. / c, 0., -a * c2;

        const double val = (a - c) / c;
        const double g2 = Math<double>::mollifier_grad(val / eps) / eps;
        const double h2 = Math<double>::mollifier_hess(val / eps) / eps / eps;

        return std::make_tuple(
            Math<double>::mollifier(val / eps), g1 * g2,
            h1 * g2 + g1 * h2 * g1.transpose());
    }
} // namespace

/// @brief Compute the gradient of the mollifier function wrt. 4 edge points and the distance squared
GradType<13> edge_edge_mollifier_gradient(
    const Eigen::Ref<const Vector3<double>>& ea0,
    const Eigen::Ref<const Vector3<double>>& ea1,
    const Eigen::Ref<const Vector3<double>>& eb0,
    const Eigen::Ref<const Vector3<double>>& eb1,
    const std::array<HEAVISIDE_TYPE, 4>& mtypes,
    const double& dist_sqr)
{
    Vector<double, 13> input;
    input << ea0, ea1, eb0, eb1, dist_sqr;

    Eigen::Matrix<int, 4, 3> vert_indices;
    vert_indices << 0, 2, 3, 1, 2, 3, 2, 0, 1, 3, 0, 1;

    Eigen::Matrix<int, 4, 9> dof_indices;
    for (int i = 0; i < 4; i++)
        dof_indices.row(i) =
            get_indices(
                vert_indices(i, 0), vert_indices(i, 1), vert_indices(i, 2))
                .transpose();

    // derivatives of point edge distances
    Eigen::Vector4d point_edge_dists;
    Eigen::Matrix<double, 4, 12> point_edge_dists_grad;
    point_edge_dists_grad.setZero();
    for (int i = 0; i < 4; i++) {
        point_edge_dists(i) = point_edge_distance(
            input.segment<3>(3 * vert_indices(i, 0)),
            input.segment<3>(3 * vert_indices(i, 1)),
            input.segment<3>(3 * vert_indices(i, 2)));
        point_edge_dists_grad(i, dof_indices.row(i)) =
            point_edge_distance_gradient(
                input.segment<3>(3 * vert_indices(i, 0)),
                input.segment<3>(3 * vert_indices(i, 1)),
                input.segment<3>(3 * vert_indices(i, 2)));
    }

    // derivatives of edge lengths
    Vector2d edge_lengths;
    Eigen::Matrix<double, 2, 12> edge_lengths_grad;
    edge_lengths_grad.setZero();
    for (int i = 0; i < 2; i++) {
        Vector3d edge = input.segment<3>(i * 6) - input.segment<3>(i * 6 + 3);
        edge_lengths(i) = edge.squaredNorm();
        edge_lengths_grad.block<1, 3>(i, 0 + 6 * i) = 2 * edge;
        edge_lengths_grad.block<1, 3>(i, 3 + 6 * i) = -2 * edge;
    }

    // derivatives of mollifier and ratios, input order : [dist_sqr,
    // edge_lengths (1 x 2), point_edge_dists (1 x 4)]
    Eigen::Vector4d mollifier;
    Eigen::Matrix<double, 4, 7> mollifier_grad;
    mollifier_grad.setZero();
    for (int i = 0; i < 4; i++) {
        const double len = edge_lengths(1 - i / 2);
        const auto [val2, g2] = func_aux_grad(
            point_edge_dists(i), len, input(12), mollifier_threshold_eps);

        Eigen::Vector3i ind;
        ind << 3 + i, 2 - i / 2, 0;
        mollifier(i) = val2;
        mollifier_grad(i, ind) = g2;
    }

    // partial products of mollifiers
    Vector4d partial_products_1;
    partial_products_1.setOnes();
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i != j)
                partial_products_1(j) *= mollifier(i);
        }
    }

    // derivatives of mollifier products, input order : [dist_sqr, edge_lengths
    // (1 x 2), point_edge_dists (1 x 4)]
    Vector<double, 7> product_grad =
        mollifier_grad.transpose() * partial_products_1;

    // derivatives wrt. input : [ea0, ea1, eb0, eb1, dist_sqr]
    Vector<double, 13> grad;
    grad << edge_lengths_grad.transpose() * product_grad.segment<2>(1)
            + point_edge_dists_grad.transpose() * product_grad.segment<4>(3),
        product_grad(0);

    return std::make_tuple(mollifier.prod(), grad);
}

/// @brief Compute the hessian of the mollifier function wrt. 4 edge points and the distance squared
HessianType<13> edge_edge_mollifier_hessian(
    const Eigen::Ref<const Vector3<double>>& ea0,
    const Eigen::Ref<const Vector3<double>>& ea1,
    const Eigen::Ref<const Vector3<double>>& eb0,
    const Eigen::Ref<const Vector3<double>>& eb1,
    const std::array<HEAVISIDE_TYPE, 4>& mtypes,
    const double& dist_sqr)
{
    Vector<double, 13> input;
    input << ea0, ea1, eb0, eb1, dist_sqr;

    Eigen::Matrix<int, 4, 3> vert_indices;
    vert_indices << 0, 2, 3, 1, 2, 3, 2, 0, 1, 3, 0, 1;

    Eigen::Matrix<int, 4, 9> dof_indices;
    for (int i = 0; i < 4; i++)
        dof_indices.row(i) =
            get_indices(
                vert_indices(i, 0), vert_indices(i, 1), vert_indices(i, 2))
                .transpose();

    // derivatives of point edge distances
    Eigen::Vector4d point_edge_dists;
    Eigen::Matrix<double, 4, 12> point_edge_dists_grad;
    point_edge_dists_grad.setZero();
    std::array<Eigen::Matrix<double, 12, 12>, 4> point_edge_dists_hess;
    for (auto& mat : point_edge_dists_hess)
        mat.setZero();
    for (int i = 0; i < 4; i++) {
        point_edge_dists(i) = point_edge_distance(
            input.segment<3>(3 * vert_indices(i, 0)),
            input.segment<3>(3 * vert_indices(i, 1)),
            input.segment<3>(3 * vert_indices(i, 2)));
        point_edge_dists_grad(i, dof_indices.row(i)) =
            point_edge_distance_gradient(
                input.segment<3>(3 * vert_indices(i, 0)),
                input.segment<3>(3 * vert_indices(i, 1)),
                input.segment<3>(3 * vert_indices(i, 2)));
        point_edge_dists_hess[i](dof_indices.row(i), dof_indices.row(i)) =
            point_edge_distance_hessian(
                input.segment<3>(3 * vert_indices(i, 0)),
                input.segment<3>(3 * vert_indices(i, 1)),
                input.segment<3>(3 * vert_indices(i, 2)));
    }

    // derivatives of edge lengths
    Vector2d edge_lengths;
    Eigen::Matrix<double, 2, 12> edge_lengths_grad;
    edge_lengths_grad.setZero();
    for (int i = 0; i < 2; i++) {
        const Vector3d edge =
            input.segment<3>(i * 6) - input.segment<3>(i * 6 + 3);
        edge_lengths(i) = edge.squaredNorm();
        edge_lengths_grad.block<1, 3>(i, 0 + 6 * i) = 2 * edge;
        edge_lengths_grad.block<1, 3>(i, 3 + 6 * i) = -2 * edge;
    }

    // derivatives of mollifier and ratios, input order : [dist_sqr,
    // edge_lengths (1 x 2), point_edge_dists (1 x 4)]
    Eigen::Vector4d mollifier;
    Eigen::Matrix<double, 4, 7> mollifier_grad;
    mollifier_grad.setZero();
    std::array<Eigen::Matrix<double, 7, 7>, 4> mollifier_hess;
    for (auto& mat : mollifier_hess)
        mat.setZero();
    for (int i = 0; i < 4; i++) {
        const double len = edge_lengths(1 - i / 2);
        const auto [val, g, h] = func_aux_hess(
            point_edge_dists(i), len, input(12), mollifier_threshold_eps);

        Eigen::Vector3i ind;
        ind << 3 + i, 2 - i / 2, 0;
        mollifier(i) = val;
        mollifier_grad(i, ind) = g;
        mollifier_hess[i](ind, ind) = h;
    }

    // partial products of mollifiers
    Vector4d partial_products_1;
    partial_products_1.setOnes();
    Eigen::Matrix<double, 4, 4> partial_products_2;
    partial_products_2.setOnes();
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i != j) {
                partial_products_1(j) *= mollifier(i);

                for (int k = 0; k < 4; k++)
                    if (i != k)
                        partial_products_2(j, k) *= mollifier(i);
            }
        }
    }
    partial_products_2.diagonal().setZero();

    // derivatives of mollifier products, input order : [dist_sqr, edge_lengths
    // (1 x 2), point_edge_dists (1 x 4)]
    Vector<double, 7> product_grad =
        mollifier_grad.transpose() * partial_products_1;
    Eigen::Matrix<double, 7, 7> product_hess =
        mollifier_grad.transpose() * partial_products_2 * mollifier_grad;
    for (int i = 0; i < 4; i++)
        product_hess += mollifier_hess[i] * partial_products_1(i);

    // derivatives wrt. input : [ea0, ea1, eb0, eb1, dist_sqr]
    Vector<double, 13> grad = Vector<double, 13>::Zero();
    Eigen::Matrix<double, 13, 13> hess = Eigen::Matrix<double, 13, 13>::Zero();
    {
        Eigen::Matrix<double, 6, 12> grads;
        grads << edge_lengths_grad, point_edge_dists_grad;

        grad(12) = product_grad(0);
        grad.head<12>() = grads.transpose() * product_grad.segment<6>(1);

        hess(12, 12) = product_hess(0, 0);
        hess.block<1, 12>(12, 0) += product_hess.block<1, 6>(0, 1) * grads;
        hess.block<12, 1>(0, 12) +=
            grads.transpose() * product_hess.block<6, 1>(1, 0);

        for (int j = 0; j < 2; j++) {
            const double val = 2 * product_grad(1 + j);
            hess.block<6, 6>(j * 6, j * 6).diagonal().array() += val;
            hess.block<3, 3>(j * 6, j * 6 + 3).diagonal().array() -= val;
            hess.block<3, 3>(j * 6 + 3, j * 6).diagonal().array() -= val;
        }
        for (int i = 0; i < 4; i++)
            hess.topLeftCorner<12, 12>() +=
                product_grad(3 + i) * point_edge_dists_hess[i];

        hess.topLeftCorner<12, 12>() +=
            grads.transpose() * product_hess.block<6, 6>(1, 1) * grads;
    }
    return std::make_tuple(mollifier.prod(), grad, hess);
}

std::array<HEAVISIDE_TYPE, 4> edge_edge_mollifier_type(
    const Eigen::Ref<const Vector3<double>>& ea0,
    const Eigen::Ref<const Vector3<double>>& ea1,
    const Eigen::Ref<const Vector3<double>>& eb0,
    const Eigen::Ref<const Vector3<double>>& eb1,
    const double& dist_sqr)
{
    std::array<HEAVISIDE_TYPE, 4> mtypes;
    mtypes[0] = (point_edge_distance(ea0, eb0, eb1) - dist_sqr)
            >= dist_sqr * mollifier_threshold_eps
        ? HEAVISIDE_TYPE::ONE
        : HEAVISIDE_TYPE::VARIANT;
    mtypes[1] = (point_edge_distance(ea1, eb0, eb1) - dist_sqr)
            >= dist_sqr * mollifier_threshold_eps
        ? HEAVISIDE_TYPE::ONE
        : HEAVISIDE_TYPE::VARIANT;
    mtypes[2] = (point_edge_distance(eb0, ea0, ea1) - dist_sqr)
            >= dist_sqr * mollifier_threshold_eps
        ? HEAVISIDE_TYPE::ONE
        : HEAVISIDE_TYPE::VARIANT;
    mtypes[3] = (point_edge_distance(eb1, ea0, ea1) - dist_sqr)
            >= dist_sqr * mollifier_threshold_eps
        ? HEAVISIDE_TYPE::ONE
        : HEAVISIDE_TYPE::VARIANT;
    return mtypes;
}

/// @brief Compute the gradient of the mollifier function wrt. 4 edge points and the distance squared
GradType<13> point_face_mollifier_gradient(
    const Eigen::Ref<const Vector3d>& p,
    const Eigen::Ref<const Vector3d>& e0,
    const Eigen::Ref<const Vector3d>& e1,
    const Eigen::Ref<const Vector3d>& e2,
    const double& dist_sqr)
{
    Vector12d x;
    x << p, e0, e1, e2;

    Vector3d vals;
    Eigen::Matrix<double, 3, 13> grads;
    for (int i = 0; i < 3; i++) {
        const int ei = i * 3 + 3;
        const int ej = ((i + 1) % 3) * 3 + 3;
        const int ek = ((i + 2) % 3) * 3 + 3;

        Vector<int, 9> ind;

        Eigen::Matrix<double, 2, 13> dist_grad =
            Eigen::Matrix<double, 2, 13>::Zero();

        const double point_edge_dist =
            point_line_distance(p, x.segment<3>(ei), x.segment<3>(ej));
        ind << 0, 1, 2, ei, ei + 1, ei + 2, ej, ej + 1, ej + 2;
        dist_grad(0, ind) =
            point_line_distance_gradient(p, x.segment<3>(ei), x.segment<3>(ej));

        const double vert_edge_dist = point_line_distance(
            x.segment<3>(ek), x.segment<3>(ei), x.segment<3>(ej));
        ind << ek, ek + 1, ek + 2, ei, ei + 1, ei + 2, ej, ej + 1, ej + 2;
        dist_grad(1, ind) = point_line_distance_gradient(
            x.segment<3>(ek), x.segment<3>(ei), x.segment<3>(ej));

        Vector3d tmp_grad;
        std::tie(vals(i), tmp_grad) = func_aux_grad(
            point_edge_dist, vert_edge_dist, dist_sqr, mollifier_threshold_eps);

        grads.row(i) = tmp_grad.head<2>().transpose() * dist_grad;
        grads(i, 12) += tmp_grad(2);
    }

    Vector<double, 13> grad = (vals(0) * vals(1)) * grads.row(2)
        + (vals(0) * vals(2)) * grads.row(1)
        + (vals(1) * vals(2)) * grads.row(0);

    return std::make_tuple(vals.prod(), grad);
}

/// @brief Compute the hessian of the mollifier function wrt. 4 edge points and the distance squared
HessianType<13> point_face_mollifier_hessian(
    const Eigen::Ref<const Vector3d>& p,
    const Eigen::Ref<const Vector3d>& e0,
    const Eigen::Ref<const Vector3d>& e1,
    const Eigen::Ref<const Vector3d>& e2,
    const double& dist_sqr)
{
    Vector12d x;
    x << p, e0, e1, e2;

    Vector3d vals;
    Eigen::Matrix<double, 3, 13> grads;
    std::array<Eigen::Matrix<double, 13, 13>, 3> hesses;
    for (int i = 0; i < 3; i++) {
        const int ei = i * 3 + 3;
        const int ej = ((i + 1) % 3) * 3 + 3;
        const int ek = ((i + 2) % 3) * 3 + 3;

        Vector<int, 9> ind;

        Eigen::Matrix<double, 2, 13> dist_grad =
            Eigen::Matrix<double, 2, 13>::Zero();
        Eigen::Matrix<double, 13, 13> point_edge_dist_hess =
                                          Eigen::Matrix<double, 13, 13>::Zero(),
                                      vert_edge_dist_hess =
                                          Eigen::Matrix<double, 13, 13>::Zero();

        const double point_edge_dist =
            point_line_distance(p, x.segment<3>(ei), x.segment<3>(ej));
        ind << 0, 1, 2, ei, ei + 1, ei + 2, ej, ej + 1, ej + 2;
        dist_grad(0, ind) =
            point_line_distance_gradient(p, x.segment<3>(ei), x.segment<3>(ej));
        point_edge_dist_hess(ind, ind) =
            point_line_distance_hessian(p, x.segment<3>(ei), x.segment<3>(ej));

        const double vert_edge_dist = point_line_distance(
            x.segment<3>(ek), x.segment<3>(ei), x.segment<3>(ej));
        ind << ek, ek + 1, ek + 2, ei, ei + 1, ei + 2, ej, ej + 1, ej + 2;
        dist_grad(1, ind) = point_line_distance_gradient(
            x.segment<3>(ek), x.segment<3>(ei), x.segment<3>(ej));
        vert_edge_dist_hess(ind, ind) = point_line_distance_hessian(
            x.segment<3>(ek), x.segment<3>(ei), x.segment<3>(ej));

        Vector3d tmp_grad;
        Matrix3d tmp_hess;
        std::tie(vals(i), tmp_grad, tmp_hess) = func_aux_hess(
            point_edge_dist, vert_edge_dist, dist_sqr, mollifier_threshold_eps);

        grads.row(i) = tmp_grad.head<2>().transpose() * dist_grad;
        grads(i, 12) += tmp_grad(2);

        hesses[i] = tmp_grad(0) * point_edge_dist_hess
            + tmp_grad(1) * vert_edge_dist_hess;
        hesses[i] +=
            dist_grad.transpose() * tmp_hess.topLeftCorner<2, 2>() * dist_grad;
        hesses[i](12, 12) += tmp_hess(2, 2);
        hesses[i].col(12) += dist_grad.transpose() * tmp_hess.block<2, 1>(0, 2);
        hesses[i].row(12) += tmp_hess.block<1, 2>(2, 0) * dist_grad;
    }

    Vector<double, 13> grad = (vals(0) * vals(1)) * grads.row(2)
        + (vals(0) * vals(2)) * grads.row(1)
        + (vals(1) * vals(2)) * grads.row(0);
    Eigen::Matrix<double, 13, 13> hess = (vals(0) * vals(1)) * hesses[2]
        + (vals(0) * vals(2)) * hesses[1] + (vals(1) * vals(2)) * hesses[0]
        + grads.row(0).transpose() * vals(2) * grads.row(1)
        + grads.row(1).transpose() * vals(2) * grads.row(0)
        + grads.row(0).transpose() * vals(1) * grads.row(2)
        + grads.row(2).transpose() * vals(1) * grads.row(0)
        + grads.row(1).transpose() * vals(0) * grads.row(2)
        + grads.row(2).transpose() * vals(0) * grads.row(1);

    return std::make_tuple(vals.prod(), grad, hess);
}

} // namespace ipc

#endif