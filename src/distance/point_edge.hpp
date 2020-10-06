#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_point.hpp>

namespace ipc {

/// @brief Compute the distance between a point and edge in 2D or 3D.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0,e1 The points of the edge.
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
auto point_edge_distance(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
{
    return point_edge_distance(p, e0, e1, point_edge_distance_type(p, e0, e1));
}

/// @brief Compute the distance between a point and edge in 2D or 3D.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0,e1 The points of the edge.
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
auto point_edge_distance(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1,
    const PointEdgeDistanceType dtype)
{
    switch (dtype) {
    case PointEdgeDistanceType::P_E0:
        return point_point_distance(p, e0);
    case PointEdgeDistanceType::P_E1:
        return point_point_distance(p, e1);
    case PointEdgeDistanceType::P_E:
        return point_line_distance(p, e0, e1);
    }
}

template <
    typename DerivedP,
    typename DerivedE0,
    typename DerivedE1,
    typename DerivedGrad>
void point_edge_distance_gradient(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1,
    Eigen::PlainObjectBase<DerivedGrad>& grad)
{
    return point_edge_distance_gradient(
        p, e0, e1, point_edge_distance_type(p, e0, e1), grad);
}

template <
    typename DerivedP,
    typename DerivedE0,
    typename DerivedE1,
    typename DerivedGrad>
void point_edge_distance_gradient(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1,
    const PointEdgeDistanceType dtype,
    Eigen::PlainObjectBase<DerivedGrad>& grad)
{
    int dim = p.size();
    assert(e0.size() == dim);
    assert(e1.size() == dim);

    grad.resize(3 * dim);
    grad.setZero();

    Eigen::VectorXd local_grad;
    switch (dtype) {
    case PointEdgeDistanceType::P_E0:
        point_point_distance_gradient(p, e0, local_grad);
        grad.head(2 * dim) = local_grad;
        break;

    case PointEdgeDistanceType::P_E1:
        point_point_distance_gradient(p, e1, local_grad);
        grad.head(dim) = local_grad.head(dim);
        grad.tail(dim) = local_grad.tail(dim);
        break;

    case PointEdgeDistanceType::P_E:
        point_line_distance_gradient(p, e0, e1, grad);
        break;
    }
}

template <
    typename DerivedP,
    typename DerivedE0,
    typename DerivedE1,
    typename DerivedHess>
void point_edge_distance_hessian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1,
    Eigen::PlainObjectBase<DerivedHess>& hess,
    bool project_to_psd = false)
{
    return point_edge_distance_hessian(
        p, e0, e1, point_edge_distance_type(p, e0, e1), hess, project_to_psd);
}

template <
    typename DerivedP,
    typename DerivedE0,
    typename DerivedE1,
    typename DerivedHess>
void point_edge_distance_hessian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1,
    const PointEdgeDistanceType dtype,
    Eigen::PlainObjectBase<DerivedHess>& hess,
    bool project_to_psd = false)
{
    int dim = p.size();
    assert(e0.size() == dim);
    assert(e1.size() == dim);

    hess.resize(3 * dim, 3 * dim);
    hess.setZero();

    Eigen::MatrixXd local_hess;
    switch (dtype) {
    case PointEdgeDistanceType::P_E0:
        point_point_distance_hessian(p, e0, local_hess, project_to_psd);
        hess.topLeftCorner(2 * dim, 2 * dim) = local_hess;
        break;

    case PointEdgeDistanceType::P_E1:
        point_point_distance_hessian(p, e1, local_hess, project_to_psd);
        hess.topLeftCorner(dim, dim) = local_hess.topLeftCorner(dim, dim);
        hess.topRightCorner(dim, dim) = local_hess.topRightCorner(dim, dim);
        hess.bottomLeftCorner(dim, dim) = local_hess.bottomLeftCorner(dim, dim);
        hess.bottomRightCorner(dim, dim) =
            local_hess.bottomRightCorner(dim, dim);
        break;

    case PointEdgeDistanceType::P_E:
        point_line_distance_hessian(p, e0, e1, hess, project_to_psd);
        break;
    }
}

} // namespace ipc
