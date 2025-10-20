#include "edge_edge.hpp"

#include "autogen.hpp"

#include <ipc/tangent/closest_point.hpp>

#include <unsupported/Eigen/KroneckerProduct>

namespace ipc {

std::tuple<Eigen::Vector3d, Eigen::Matrix<double, 3, 12>>
line_line_closest_point_direction_gradient(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
    const Eigen::Vector3d ea = ea1 - ea0;
    const Eigen::Vector3d eb = eb1 - eb0;
    const Eigen::Vector3d normal = ea.cross(eb);
    const Eigen::Matrix<double, 3, 6> cross_grad =
        cross_product_gradient(ea, eb);

    const double normal_sqr_norm = normal.squaredNorm();
    const Eigen::Vector3d t = eb0 - ea0;
    const Eigen::Vector3d vec = (t.dot(normal) / normal_sqr_norm) * normal;

    Eigen::Matrix<double, 3, 12> grad = Eigen::Matrix<double, 3, 12>::Zero();
    // derivative wrt. normal, tempararily store here
    grad(Eigen::all, { 3, 4, 5 }) =
        (normal * t.transpose()
         + t.transpose() * normal * Eigen::Matrix3d::Identity()
         - 2 * vec * normal.transpose())
        / normal_sqr_norm;
    // derivative wrt. t
    grad.middleCols<3>(6) = normal * normal.transpose() / normal_sqr_norm;
    grad.leftCols<3>() = -grad.middleCols<3>(6);
    // contributions from n
    grad(Eigen::all, { 3, 4, 5, 9, 10, 11 }) =
        grad(Eigen::all, { 3, 4, 5 }) * cross_grad;
    grad(Eigen::all, { 0, 1, 2, 6, 7, 8 }) -=
        grad(Eigen::all, { 3, 4, 5, 9, 10, 11 });

    return std::make_tuple(vec, grad);
}

std::tuple<
    Eigen::Vector3d,
    Eigen::Matrix<double, 3, 12>,
    std::array<Matrix12d, 3>>
line_line_closest_point_direction_hessian(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
    const Eigen::Vector3d ea = ea1 - ea0;
    const Eigen::Vector3d eb = eb1 - eb0;
    const Eigen::Vector3d normal = ea.cross(eb);
    const Eigen::Matrix<double, 3, 6> cross_grad =
        cross_product_gradient(ea, eb);
    const std::array<Matrix6d, 3> cross_hess = cross_product_hessian(ea, eb);

    const double normal_sqr_norm = normal.squaredNorm();
    const Eigen::Vector3d t = eb0 - ea0;
    const Eigen::Vector3d vec = (t.dot(normal) / normal_sqr_norm) * normal;

    Eigen::Matrix<double, 3, 12> grad = Eigen::Matrix<double, 3, 12>::Zero();
    // derivative wrt. normal and t
    Eigen::Matrix3d grad_normal =
        (normal * t.transpose() - 2 * vec * normal.transpose())
        / normal_sqr_norm;
    grad_normal.diagonal().array() +=
        (t.transpose() * normal)(0) / normal_sqr_norm;
    Eigen::Matrix3d grad_t = normal * normal.transpose() / normal_sqr_norm;
    grad.middleCols<3>(6) = grad_t;
    grad.leftCols<3>() = -grad_t;
    // contributions from n to ea, eb
    grad(Eigen::all, { 3, 4, 5, 9, 10, 11 }) = grad_normal * cross_grad;
    grad(Eigen::all, { 0, 1, 2, 6, 7, 8 }) -=
        grad(Eigen::all, { 3, 4, 5, 9, 10, 11 });

    // hessian of output wrt. normal
    std::array<Eigen::Matrix3d, 3> hess_wrt_normal;
    // hessian of output wrt. normal (row) and t (col)
    std::array<Eigen::Matrix3d, 3> mixed_hess;
    for (int d = 0; d < 3; d++) {
        hess_wrt_normal[d] = -2
            * (normal * grad_normal.row(d)
               + grad_normal.row(d).transpose() * normal.transpose());
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
    inner_grad << -cross_grad.leftCols<3>(), cross_grad.leftCols<3>(),
        -cross_grad.rightCols<3>(), cross_grad.rightCols<3>();

    std::array<Matrix12d, 3> inner_hess;
    {
        Eigen::Matrix2d tmp;
        tmp << 1, -1, -1, 1;
        for (int d = 0; d < 3; d++) {
            inner_hess[d]
                << Eigen::KroneckerProduct<Eigen::Matrix2d, Eigen::Matrix3d>(
                       tmp, cross_hess[d].block<3, 3>(0, 0)),
                Eigen::KroneckerProduct<Eigen::Matrix2d, Eigen::Matrix3d>(
                    tmp, cross_hess[d].block<3, 3>(0, 3)),
                Eigen::KroneckerProduct<Eigen::Matrix2d, Eigen::Matrix3d>(
                    tmp, cross_hess[d].block<3, 3>(3, 0)),
                Eigen::KroneckerProduct<Eigen::Matrix2d, Eigen::Matrix3d>(
                    tmp, cross_hess[d].block<3, 3>(3, 3));
        }
    }

    std::array<Matrix12d, 3> hess;
    for (int d = 0; d < 3; d++) {
        Eigen::Matrix<double, 12, 3> tmp =
            inner_grad.transpose() * mixed_hess[d];
        hess[d] = inner_grad.transpose() * hess_wrt_normal[d] * inner_grad;
        hess[d].middleCols<3>(6) += tmp;
        hess[d].middleCols<3>(0) -= tmp;
        hess[d].middleRows<3>(6) += tmp.transpose();
        hess[d].middleRows<3>(0) -= tmp.transpose();
        for (int i = 0; i < 3; i++) {
            hess[d] += inner_hess[i] * grad_normal(d, i);
        }
    }

    return std::make_tuple(vec, grad, hess);
}

std::tuple<Vector6d, Eigen::Matrix<double, 6, 12>>
line_line_closest_point_pairs_gradient(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
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

    return { out, grad };
}

std::tuple<Vector6d, Eigen::Matrix<double, 6, 12>, std::array<Matrix12d, 6>>
line_line_closest_point_pairs_hessian(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
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
    for (auto& h : hess) {
        h.setZero();
    }

    {
        Eigen::Matrix<double, 12, 12> Ha, Hb;
        autogen::edge_edge_closest_point_hessian_a(
            ea0.x(), ea0.y(), ea0.z(), ea1.x(), ea1.y(), ea1.z(), eb0.x(),
            eb0.y(), eb0.z(), eb1.x(), eb1.y(), eb1.z(), Ha.data());
        autogen::edge_edge_closest_point_hessian_b(
            ea0.x(), ea0.y(), ea0.z(), ea1.x(), ea1.y(), ea1.z(), eb0.x(),
            eb0.y(), eb0.z(), eb1.x(), eb1.y(), eb1.z(), Hb.data());

        for (int d = 0; d < 3; d++) {
            // wrt. uv
            hess[d] += (ea1(d) - ea0(d)) * Ha;
            hess[d + 3] += (eb1(d) - eb0(d)) * Hb;

            // wrt. ea0 & uv(0)
            hess[d].row(d) -= J.row(0);
            hess[d].col(d) -= J.row(0).transpose();

            // wrt. ea1 & uv(0)
            hess[d].row(d + 3) += J.row(0);
            hess[d].col(d + 3) += J.row(0).transpose();

            // wrt. eb0 & uv(1)
            hess[d + 3].row(d + 6) -= J.row(1);
            hess[d + 3].col(d + 6) -= J.row(1).transpose();

            // wrt. eb1 & uv(1)
            hess[d + 3].row(d + 9) += J.row(1);
            hess[d + 3].col(d + 9) += J.row(1).transpose();
        }
    }

    return { out, grad, hess };
}

template <typename scalar>
scalar line_line_sqr_distance(
    Eigen::ConstRef<Eigen::Vector3<scalar>> ea0,
    Eigen::ConstRef<Eigen::Vector3<scalar>> ea1,
    Eigen::ConstRef<Eigen::Vector3<scalar>> eb0,
    Eigen::ConstRef<Eigen::Vector3<scalar>> eb1)
{
    const Eigen::Vector3<scalar> normal = (ea1 - ea0).cross(eb1 - eb0);
    const scalar line_to_line = (eb0 - ea0).dot(normal);
    return line_to_line * line_to_line / normal.squaredNorm();
}

template <typename scalar>
scalar edge_edge_sqr_distance(
    Eigen::ConstRef<Eigen::Vector3<scalar>> ea0,
    Eigen::ConstRef<Eigen::Vector3<scalar>> ea1,
    Eigen::ConstRef<Eigen::Vector3<scalar>> eb0,
    Eigen::ConstRef<Eigen::Vector3<scalar>> eb1,
    EdgeEdgeDistanceType dtype)
{
    if constexpr (std::is_same<double, scalar>::value) {
        if (dtype == EdgeEdgeDistanceType::AUTO) {
            dtype = edge_edge_distance_type(ea0, ea1, eb0, eb1);
        }
    }

    switch (dtype) {
    case EdgeEdgeDistanceType::EA0_EB0:
        return PointEdgeDistance<scalar, 3>::point_point_sqr_distance(ea0, eb0);

    case EdgeEdgeDistanceType::EA0_EB1:
        return PointEdgeDistance<scalar, 3>::point_point_sqr_distance(ea0, eb1);

    case EdgeEdgeDistanceType::EA1_EB0:
        return PointEdgeDistance<scalar, 3>::point_point_sqr_distance(ea1, eb0);

    case EdgeEdgeDistanceType::EA1_EB1:
        return PointEdgeDistance<scalar, 3>::point_point_sqr_distance(ea1, eb1);

    case EdgeEdgeDistanceType::EA_EB0:
        return PointEdgeDistance<scalar, 3>::point_line_sqr_distance(
            eb0, ea0, ea1);

    case EdgeEdgeDistanceType::EA_EB1:
        return PointEdgeDistance<scalar, 3>::point_line_sqr_distance(
            eb1, ea0, ea1);

    case EdgeEdgeDistanceType::EA0_EB:
        return PointEdgeDistance<scalar, 3>::point_line_sqr_distance(
            ea0, eb0, eb1);

    case EdgeEdgeDistanceType::EA1_EB:
        return PointEdgeDistance<scalar, 3>::point_line_sqr_distance(
            ea1, eb0, eb1);

    case EdgeEdgeDistanceType::EA_EB:
        return line_line_sqr_distance<scalar>(ea0, ea1, eb0, eb1);

    default:
        throw std::invalid_argument(
            "Invalid distance type for edge-edge distance!");
    }
}

template <typename scalar>
Eigen::Vector3<scalar> line_line_closest_point_direction(
    Eigen::ConstRef<Eigen::Vector3<scalar>> ea0,
    Eigen::ConstRef<Eigen::Vector3<scalar>> ea1,
    Eigen::ConstRef<Eigen::Vector3<scalar>> eb0,
    Eigen::ConstRef<Eigen::Vector3<scalar>> eb1)
{
    const Eigen::Vector3<scalar> normal = (ea1 - ea0).cross(eb1 - eb0);
    return ((eb0 - ea0).dot(normal) / normal.squaredNorm()) * normal;
}

template <typename scalar>
Eigen::Matrix<scalar, 3, 2> line_line_closest_point_pairs(
    Eigen::ConstRef<Eigen::Vector3<scalar>> ea0,
    Eigen::ConstRef<Eigen::Vector3<scalar>> ea1,
    Eigen::ConstRef<Eigen::Vector3<scalar>> eb0,
    Eigen::ConstRef<Eigen::Vector3<scalar>> eb1)
{
    const Eigen::Vector3<scalar> ta = ea1 - ea0;
    const Eigen::Vector3<scalar> tb = eb1 - eb0;
    const scalar la = ta.squaredNorm();
    const scalar lb = tb.squaredNorm();
    const scalar lab = ta.dot(tb);
    const Eigen::Vector3<scalar> d = eb0 - ea0;

    Eigen::Matrix<scalar, 3, 2> out;
    const scalar fac = la * lb - pow(lab, 2);
    out.col(0) = ea0 + (lb * ta.dot(d) - lab * tb.dot(d)) / fac * ta;
    out.col(1) = eb0 + (lab * ta.dot(d) - la * tb.dot(d)) / fac * tb;

    return out;
}

/// @brief Computes the direction of the closest point pair
/// @param ea0 Vertex 0 of edge 0
/// @param ea1 Vertex 1 of edge 0
/// @param eb0 Vertex 0 of edge 1
/// @param eb1 Vertex 1 of edge 1
/// @param dtype Edge-edge distance type
/// @return Difference of the pair of closest point, pointing from edge 0 to edge 1
template <typename scalar>
Eigen::Vector3<scalar> edge_edge_closest_point_direction(
    Eigen::ConstRef<Eigen::Vector3<scalar>> ea0,
    Eigen::ConstRef<Eigen::Vector3<scalar>> ea1,
    Eigen::ConstRef<Eigen::Vector3<scalar>> eb0,
    Eigen::ConstRef<Eigen::Vector3<scalar>> eb1,
    EdgeEdgeDistanceType dtype)
{
    if constexpr (std::is_same<double, scalar>::value) {
        if (dtype == EdgeEdgeDistanceType::AUTO) {
            dtype = edge_edge_distance_type(ea0, ea1, eb0, eb1);
        }
    }

    switch (dtype) {
    case EdgeEdgeDistanceType::EA0_EB0:
        return (eb0 - ea0);

    case EdgeEdgeDistanceType::EA0_EB1:
        return (eb1 - ea0);

    case EdgeEdgeDistanceType::EA1_EB0:
        return (eb0 - ea1);

    case EdgeEdgeDistanceType::EA1_EB1:
        return (eb1 - ea1);

    case EdgeEdgeDistanceType::EA_EB0:
        return PointEdgeDistance<scalar, 3>::point_line_closest_point_direction(
            eb0, ea0, ea1);

    case EdgeEdgeDistanceType::EA_EB1:
        return PointEdgeDistance<scalar, 3>::point_line_closest_point_direction(
            eb1, ea0, ea1);

    case EdgeEdgeDistanceType::EA0_EB:
        return -PointEdgeDistance<
            scalar, 3>::point_line_closest_point_direction(ea0, eb0, eb1);

    case EdgeEdgeDistanceType::EA1_EB:
        return -PointEdgeDistance<
            scalar, 3>::point_line_closest_point_direction(ea1, eb0, eb1);

    case EdgeEdgeDistanceType::EA_EB:
        return line_line_closest_point_direction<scalar>(ea0, ea1, eb0, eb1);

    default:
        throw std::invalid_argument(
            "Invalid distance type for edge-edge distance!");
    }
}

template <typename scalar>
Eigen::Matrix<scalar, 3, 2> edge_edge_closest_point_pairs(
    Eigen::ConstRef<Eigen::Vector3<scalar>> ea0,
    Eigen::ConstRef<Eigen::Vector3<scalar>> ea1,
    Eigen::ConstRef<Eigen::Vector3<scalar>> eb0,
    Eigen::ConstRef<Eigen::Vector3<scalar>> eb1,
    EdgeEdgeDistanceType dtype)
{
    if constexpr (std::is_same<double, scalar>::value) {
        if (dtype == EdgeEdgeDistanceType::AUTO) {
            dtype = edge_edge_distance_type(ea0, ea1, eb0, eb1);
        }
    }

    Eigen::Matrix<scalar, 3, 2> out;
    switch (dtype) {
    case EdgeEdgeDistanceType::EA0_EB0:
        out << ea0, eb0;
        break;

    case EdgeEdgeDistanceType::EA0_EB1:
        out << ea0, eb1;
        break;

    case EdgeEdgeDistanceType::EA1_EB0:
        out << ea1, eb0;
        break;

    case EdgeEdgeDistanceType::EA1_EB1:
        out << ea1, eb1;
        break;

    case EdgeEdgeDistanceType::EA_EB0:
        out << eb0
                - PointEdgeDistance<scalar, 3>::
                    point_line_closest_point_direction(eb0, ea0, ea1),
            eb0;
        break;

    case EdgeEdgeDistanceType::EA_EB1:
        out << eb1
                - PointEdgeDistance<scalar, 3>::
                    point_line_closest_point_direction(eb1, ea0, ea1),
            eb1;
        break;

    case EdgeEdgeDistanceType::EA0_EB:
        out << ea0,
            ea0
            - PointEdgeDistance<scalar, 3>::point_line_closest_point_direction(
                ea0, eb0, eb1);
        break;

    case EdgeEdgeDistanceType::EA1_EB:
        out << ea1,
            ea1
            - PointEdgeDistance<scalar, 3>::point_line_closest_point_direction(
                ea1, eb0, eb1);
        break;

    case EdgeEdgeDistanceType::EA_EB:
        out = line_line_closest_point_pairs<scalar>(ea0, ea1, eb0, eb1);
        break;

    default:
        throw std::invalid_argument(
            "Invalid distance type for edge-edge distance!");
    }

    return out;
}

template Eigen::Vector3d edge_edge_closest_point_direction(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1,
    EdgeEdgeDistanceType dtype);
template Eigen::Vector3<ADGrad<9>> edge_edge_closest_point_direction(
    Eigen::ConstRef<Eigen::Vector3<ADGrad<9>>> ea0,
    Eigen::ConstRef<Eigen::Vector3<ADGrad<9>>> ea1,
    Eigen::ConstRef<Eigen::Vector3<ADGrad<9>>> eb0,
    Eigen::ConstRef<Eigen::Vector3<ADGrad<9>>> eb1,
    EdgeEdgeDistanceType dtype);
template Eigen::Vector3<ADHessian<9>> edge_edge_closest_point_direction(
    Eigen::ConstRef<Eigen::Vector3<ADHessian<9>>> ea0,
    Eigen::ConstRef<Eigen::Vector3<ADHessian<9>>> ea1,
    Eigen::ConstRef<Eigen::Vector3<ADHessian<9>>> eb0,
    Eigen::ConstRef<Eigen::Vector3<ADHessian<9>>> eb1,
    EdgeEdgeDistanceType dtype);
template Eigen::Vector3<ADGrad<12>> edge_edge_closest_point_direction(
    Eigen::ConstRef<Eigen::Vector3<ADGrad<12>>> ea0,
    Eigen::ConstRef<Eigen::Vector3<ADGrad<12>>> ea1,
    Eigen::ConstRef<Eigen::Vector3<ADGrad<12>>> eb0,
    Eigen::ConstRef<Eigen::Vector3<ADGrad<12>>> eb1,
    EdgeEdgeDistanceType dtype);
template Eigen::Vector3<ADHessian<12>> edge_edge_closest_point_direction(
    Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>> ea0,
    Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>> ea1,
    Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>> eb0,
    Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>> eb1,
    EdgeEdgeDistanceType dtype);
template Eigen::Vector3<ADGrad<13>> edge_edge_closest_point_direction(
    Eigen::ConstRef<Eigen::Vector3<ADGrad<13>>> ea0,
    Eigen::ConstRef<Eigen::Vector3<ADGrad<13>>> ea1,
    Eigen::ConstRef<Eigen::Vector3<ADGrad<13>>> eb0,
    Eigen::ConstRef<Eigen::Vector3<ADGrad<13>>> eb1,
    EdgeEdgeDistanceType dtype);
template Eigen::Vector3<ADHessian<13>> edge_edge_closest_point_direction(
    Eigen::ConstRef<Eigen::Vector3<ADHessian<13>>> ea0,
    Eigen::ConstRef<Eigen::Vector3<ADHessian<13>>> ea1,
    Eigen::ConstRef<Eigen::Vector3<ADHessian<13>>> eb0,
    Eigen::ConstRef<Eigen::Vector3<ADHessian<13>>> eb1,
    EdgeEdgeDistanceType dtype);

template double edge_edge_sqr_distance(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1,
    EdgeEdgeDistanceType dtype);
template ADGrad<9> edge_edge_sqr_distance(
    Eigen::ConstRef<Eigen::Vector3<ADGrad<9>>> ea0,
    Eigen::ConstRef<Eigen::Vector3<ADGrad<9>>> ea1,
    Eigen::ConstRef<Eigen::Vector3<ADGrad<9>>> eb0,
    Eigen::ConstRef<Eigen::Vector3<ADGrad<9>>> eb1,
    EdgeEdgeDistanceType dtype);
template ADHessian<9> edge_edge_sqr_distance(
    Eigen::ConstRef<Eigen::Vector3<ADHessian<9>>> ea0,
    Eigen::ConstRef<Eigen::Vector3<ADHessian<9>>> ea1,
    Eigen::ConstRef<Eigen::Vector3<ADHessian<9>>> eb0,
    Eigen::ConstRef<Eigen::Vector3<ADHessian<9>>> eb1,
    EdgeEdgeDistanceType dtype);
template ADGrad<12> edge_edge_sqr_distance(
    Eigen::ConstRef<Eigen::Vector3<ADGrad<12>>> ea0,
    Eigen::ConstRef<Eigen::Vector3<ADGrad<12>>> ea1,
    Eigen::ConstRef<Eigen::Vector3<ADGrad<12>>> eb0,
    Eigen::ConstRef<Eigen::Vector3<ADGrad<12>>> eb1,
    EdgeEdgeDistanceType dtype);
template ADHessian<12> edge_edge_sqr_distance(
    Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>> ea0,
    Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>> ea1,
    Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>> eb0,
    Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>> eb1,
    EdgeEdgeDistanceType dtype);
template ADGrad<13> edge_edge_sqr_distance(
    Eigen::ConstRef<Eigen::Vector3<ADGrad<13>>> ea0,
    Eigen::ConstRef<Eigen::Vector3<ADGrad<13>>> ea1,
    Eigen::ConstRef<Eigen::Vector3<ADGrad<13>>> eb0,
    Eigen::ConstRef<Eigen::Vector3<ADGrad<13>>> eb1,
    EdgeEdgeDistanceType dtype);
template ADHessian<13> edge_edge_sqr_distance(
    Eigen::ConstRef<Eigen::Vector3<ADHessian<13>>> ea0,
    Eigen::ConstRef<Eigen::Vector3<ADHessian<13>>> ea1,
    Eigen::ConstRef<Eigen::Vector3<ADHessian<13>>> eb0,
    Eigen::ConstRef<Eigen::Vector3<ADHessian<13>>> eb1,
    EdgeEdgeDistanceType dtype);

template Eigen::Matrix<double, 3, 2> line_line_closest_point_pairs(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1);
template Eigen::Matrix<ADGrad<12>, 3, 2> line_line_closest_point_pairs(
    Eigen::ConstRef<Eigen::Vector3<ADGrad<12>>> ea0,
    Eigen::ConstRef<Eigen::Vector3<ADGrad<12>>> ea1,
    Eigen::ConstRef<Eigen::Vector3<ADGrad<12>>> eb0,
    Eigen::ConstRef<Eigen::Vector3<ADGrad<12>>> eb1);
template Eigen::Matrix<ADHessian<12>, 3, 2> line_line_closest_point_pairs(
    Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>> ea0,
    Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>> ea1,
    Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>> eb0,
    Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>> eb1);
} // namespace ipc
