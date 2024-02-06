#include "mollifier.hpp"
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_edge.hpp>
#include "point_edge.hpp"

namespace ipc {
namespace {
    Vector<int, 9> get_indices(int i, int j, int k)
    {
        Vector<int, 9> out;
        out << 3 * i + 0, 3 * i + 1, 3 * i + 2, 3 * j + 0, 3 * j + 1, 3 * j + 2,
            3 * k + 0, 3 * k + 1, 3 * k + 2;
        return out;
    }

    double func_aux1(const double& a, const double& b, const double& c, const double& eps)
    {
        return (a - c) / b / eps;
    }

    Eigen::Vector3d
    func_aux1_grad(const double& a, const double& b, const double& c, const double& eps)
    {
        return Eigen::Vector3d(1. / b, -func_aux1(a, b, c, 1.) / b, -1. / b) / eps;
    }

    Eigen::Matrix3d
    func_aux1_hess(const double& a, const double& b, const double& c, const double& eps)
    {
        const double b2 = 1. / b / b;
        Eigen::Matrix3d out;
        out << 0., -b2, 0., -b2, 2. * func_aux1(a, b, c, 1.) * b2, b2, 0., b2, 0.;
        return out / eps;
    }
} // namespace

/// @brief Compute the gradient of the mollifier function wrt. 4 edge points and the distance squared
std::pair<double, Vector<double, 13>> edge_edge_mollifier_grad(
    const Vector3<double>& ea0,
    const Vector3<double>& ea1,
    const Vector3<double>& eb0,
    const Vector3<double>& eb1,
    const std::array<HEAVISIDE_TYPE, 4>& mtypes,
    const double& dist_sqr)
{
    const double da0 = point_edge_distance(ea0, eb0, eb1);
    const double da1 = point_edge_distance(ea1, eb0, eb1);
    const double db0 = point_edge_distance(eb0, ea0, ea1);
    const double db1 = point_edge_distance(eb1, ea0, ea1);

    Vector12d da0_wrt_x, da1_wrt_x, db0_wrt_x, db1_wrt_x;
    {
        Vector9d tmp = point_edge_distance_gradient(ea0, eb0, eb1);
        da0_wrt_x << tmp.head(3), Eigen::Vector3d::Zero(), tmp.tail(6);

        tmp = point_edge_distance_gradient(ea1, eb0, eb1);
        da1_wrt_x << Eigen::Vector3d::Zero(), tmp;

        tmp = point_edge_distance_gradient(eb0, ea0, ea1);
        db0_wrt_x << tmp.tail(6), tmp.head(3), Eigen::Vector3d::Zero();

        tmp = point_edge_distance_gradient(eb1, ea0, ea1);
        db1_wrt_x << tmp.tail(6), Eigen::Vector3d::Zero(), tmp.head(3);
    }

    const double db = point_point_distance(eb0, eb1) * mollifier_threshold_eps;
    const double da = point_point_distance(ea0, ea1) * mollifier_threshold_eps;

    Vector6d da_wrt_x =
        point_point_distance_gradient(ea0, ea1) * mollifier_threshold_eps;
    Vector6d db_wrt_x =
        point_point_distance_gradient(eb0, eb1) * mollifier_threshold_eps;

    const double ma0 = func_aux1(da0, db, dist_sqr, 1.);
    const double ma1 = func_aux1(da1, db, dist_sqr, 1.);
    const double mb0 = func_aux1(db0, da, dist_sqr, 1.);
    const double mb1 = func_aux1(db1, da, dist_sqr, 1.);

    const double a = Math<double>::mollifier(ma0);
    const double b = Math<double>::mollifier(ma1);
    const double c = Math<double>::mollifier(mb0);
    const double d = Math<double>::mollifier(mb1);

    const double a_wrt_ma0 = Math<double>::mollifier_grad(ma0);
    const double b_wrt_ma1 = Math<double>::mollifier_grad(ma1);
    const double c_wrt_mb0 = Math<double>::mollifier_grad(mb0);
    const double d_wrt_mb1 = Math<double>::mollifier_grad(mb1);

    const double ab = a * b, cd = c * d;
    Vector12d ab_grad = a * b_wrt_ma1 * da1_wrt_x + b * a_wrt_ma0 * da0_wrt_x;
    ab_grad.tail(6) -= db_wrt_x * (a * b_wrt_ma1 * ma1 + b * a_wrt_ma0 * ma0);
    ab_grad /= db;
    Vector12d cd_grad = c * d_wrt_mb1 * db1_wrt_x + d * c_wrt_mb0 * db0_wrt_x;
    cd_grad.head(6) -= da_wrt_x * (c * d_wrt_mb1 * mb1 + d * c_wrt_mb0 * mb0);
    cd_grad /= da;

    Vector<double, 13> out;
    out.head(12) = ab_grad * cd + cd_grad * ab;
    out(12) = -(a * b_wrt_ma1 + b * a_wrt_ma0) * cd / db
        - (c * d_wrt_mb1 + d * c_wrt_mb0) * ab / da;

    return std::make_pair(ab * cd, out);
}

/// @brief Compute the hessian of the mollifier function wrt. 4 edge points and the distance squared
std::tuple<double, Vector<double, 13>, Eigen::Matrix<double, 13, 13>>
edge_edge_mollifier_hessian(
    const Vector3<double>& ea0,
    const Vector3<double>& ea1,
    const Vector3<double>& eb0,
    const Vector3<double>& eb1,
    const std::array<HEAVISIDE_TYPE, 4>& mtypes,
    const double& dist_sqr)
{
    Vector<double, 13> input;
    input << ea0, ea1, eb0, eb1, dist_sqr;

#ifndef DERIVATIVES_WITH_AUTODIFF
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
    for (int i = 0; i < 4; i++)
    {
        point_edge_dists(i) =
            point_edge_distance(input.segment<3>(3 * vert_indices(i, 0)),
            input.segment<3>(3 * vert_indices(i, 1)),
            input.segment<3>(3 * vert_indices(i, 2))); 
        point_edge_dists_grad(i, dof_indices.row(i)) =
            point_edge_distance_gradient(input.segment<3>(3 * vert_indices(i, 0)),
            input.segment<3>(3 * vert_indices(i, 1)),
            input.segment<3>(3 * vert_indices(i, 2)));
        point_edge_dists_hess[i](dof_indices.row(i), dof_indices.row(i)) =
            point_edge_distance_hessian(input.segment<3>(3 * vert_indices(i, 0)),
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

    // derivatives of mollifier and ratios, input order : [dist_sqr, edge_lengths
    // (1 x 2), point_edge_dists (1 x 4)]
    Eigen::Vector4d mollifier;
    Eigen::Matrix<double, 4, 7> mollifier_grad;
    mollifier_grad.setZero();
    std::array<Eigen::Matrix<double, 7, 7>, 4> mollifier_hess;
    for (auto& mat : mollifier_hess)
        mat.setZero();
    for (int i = 0; i < 4; i++) {
        const double len = edge_lengths(1 - i / 2);
        const double val1 = func_aux1(point_edge_dists(i), len, input(12), mollifier_threshold_eps);
        const Vector3d g1 = func_aux1_grad(point_edge_dists(i), len, input(12), mollifier_threshold_eps);
        const Eigen::Matrix3d h1 = func_aux1_hess(point_edge_dists(i), len, input(12), mollifier_threshold_eps);

        const double val2 = Math<double>::mollifier(val1);
        const double g2 = Math<double>::mollifier_grad(val1);
        const double h2 = Math<double>::mollifier_hess(val1);

        Eigen::Vector3i ind;
        ind << 3 + i, 2 - i / 2, 0;
        mollifier(i) = val2;
        mollifier_grad(i, ind) = g2 * g1;
        mollifier_hess[i](ind, ind) = g1 * h2 * g1.transpose() + h1 * g2;
    }

    // partial products of mollifiers
    Vector<double, 4> partial_products_1;
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
    Vector<double, 7> product_grad = mollifier_grad.transpose() * partial_products_1;
    Eigen::Matrix<double, 7, 7> product_hess = mollifier_grad.transpose() * partial_products_2 * mollifier_grad;
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
        hess.block<12, 1>(0, 12) += grads.transpose() * product_hess.block<6, 1>(1, 0);

        for (int j = 0; j < 2; j++)
        {
            hess.block<6, 6>(j * 6, j * 6).diagonal().array() += 2 * product_grad(1 + j);
            hess.block<3, 3>(j * 6, j * 6 + 3).diagonal().array() -= 2 * product_grad(1 + j);
            hess.block<3, 3>(j * 6 + 3, j * 6).diagonal().array() -= 2 * product_grad(1 + j);
        }
        for (int i = 0; i < 4; i++)
            hess.topLeftCorner<12, 12>() += product_grad(3 + i) * point_edge_dists_hess[i];

        hess.topLeftCorner<12, 12>() += grads.transpose() * product_hess.block<6, 6>(1, 1) * grads;
    }
    return std::make_tuple(mollifier.prod(), grad, hess);
#else
    DiffScalarBase::setVariableCount(13);
    using T = ADHessian<13>;
    Vector<T, 13> input_ad = slice_positions<T, 13, 1>(input);

    const T da = (input_ad.segment<3>(3) - input_ad.head<3>()).squaredNorm()
        * mollifier_threshold_eps;
    const T db = (input_ad.segment<3>(9) - input_ad.segment<3>(6)).squaredNorm()
        * mollifier_threshold_eps;
    T aAD = (mtypes[0] == HEAVISIDE_TYPE::VARIANT) ? Math<T>::mollifier(
                (PointEdgeDistance<T, 3>::point_edge_sqr_distance(
                     input_ad.head<3>(), input_ad.segment<3>(6),
                     input_ad.segment<3>(9))
                 - input_ad(12)) / db) : T(1.);
    T bAD = (mtypes[1] == HEAVISIDE_TYPE::VARIANT) ? Math<T>::mollifier(
                (PointEdgeDistance<T, 3>::point_edge_sqr_distance(
                     input_ad.segment<3>(3), input_ad.segment<3>(6),
                     input_ad.segment<3>(9))
                 - input_ad(12)) / db) : T(1.);
    T cAD = (mtypes[2] == HEAVISIDE_TYPE::VARIANT) ? Math<T>::mollifier(
                (PointEdgeDistance<T, 3>::point_edge_sqr_distance(
                     input_ad.segment<3>(6), input_ad.head<3>(),
                     input_ad.segment<3>(3))
                 - input_ad(12)) / da) : T(1.);
    T dAD = (mtypes[3] == HEAVISIDE_TYPE::VARIANT) ? Math<T>::mollifier(
                (PointEdgeDistance<T, 3>::point_edge_sqr_distance(
                     input_ad.segment<3>(9), input_ad.head<3>(),
                     input_ad.segment<3>(3))
                 - input_ad(12)) / da) : T(1.);

    T outAD = aAD * bAD * cAD * dAD;

    return std::make_tuple(
        outAD.getValue(), outAD.getGradient(), outAD.getHessian());
#endif
}

} // namespace ipc
