#pragma once

#include <distance/distance_type.hpp>
#include <distance/line_line.hpp>
#include <distance/point_edge.hpp>
#include <distance/point_point.hpp>

namespace ipc {

/// @brief Compute the distance between a two lines segments in 3D.
/// @note The distance is actually squared distance.
/// @param ea0,ea1 The points of the first edge.
/// @param eb0,eb1 The points of the second edge.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
auto edge_edge_distance(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
{
    switch (edge_edge_distance_type(ea0, ea1, eb0, eb1)) {
    case EdgeEdgeDistanceType::EA0_EB0:
        return point_point_distance(ea0, eb0);

    case EdgeEdgeDistanceType::EA0_EB1:
        return point_point_distance(ea0, eb1);

    case EdgeEdgeDistanceType::EA1_EB0:
        return point_point_distance(ea1, eb0);

    case EdgeEdgeDistanceType::EA1_EB1:
        return point_point_distance(ea1, eb1);

    case EdgeEdgeDistanceType::EA_EB0:
        return point_edge_distance(eb0, ea0, ea1);

    case EdgeEdgeDistanceType::EA_EB1:
        return point_edge_distance(eb1, ea0, ea1);

    case EdgeEdgeDistanceType::EA0_EB:
        return point_edge_distance(ea0, eb0, eb1);

    case EdgeEdgeDistanceType::EA1_EB:
        return point_edge_distance(ea1, eb0, eb1);

    case EdgeEdgeDistanceType::EA_EB:
        return line_line_distance(ea0, ea1, eb0, eb1);
    }

    throw "something went wrong in edge_edge_distance";
}

/// @brief Compute the gradient of the distance between a two lines segments.
/// @note The distance is actually squared distance.
/// @param[in] ea0,ea1 The points of the first edge.
/// @param[in] eb0,eb1 The points of the second edge.
/// @param[out] grad The computed gradient.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1,
    typename DerivedGrad>
auto edge_edge_distance_gradient(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1,
    Eigen::MatrixBase<DerivedGrad>& grad)
{
    int dim = ea0.size();
    assert(ea1.size() == dim);
    assert(eb0.size() == dim);
    assert(eb1.size() == dim);

    grad.resize(4 * dim);
    grad.setZero();

    Eigen::VectorXd local_grad;
    switch (edge_edge_distance_type(ea0, ea1, eb0, eb1)) {
    case EdgeEdgeDistanceType::EA0_EB0:
        point_point_distance_gradient(ea0, eb0, local_grad);
        grad.head(dim) = local_grad.head(dim);
        grad.segment(2 * dim, dim) = local_grad.tail(dim);
        break;

    case EdgeEdgeDistanceType::EA0_EB1:
        point_point_distance_gradient(ea0, eb1, local_grad);
        grad.head(dim) = local_grad.head(dim);
        grad.tail(dim) = local_grad.tail(dim);
        break;

    case EdgeEdgeDistanceType::EA1_EB0:
        point_point_distance_gradient(ea1, eb0, local_grad);
        grad.segment(dim, dim) = local_grad.head(dim);
        grad.segment(2 * dim, dim) = local_grad.tail(dim);
        break;

    case EdgeEdgeDistanceType::EA1_EB1:
        point_point_distance_gradient(ea1, eb1, local_grad);
        grad.segment(dim, dim) = local_grad.head(dim);
        grad.tail(dim) = local_grad.tail(dim);
        break;

    case EdgeEdgeDistanceType::EA_EB0:
        point_edge_distance_gradient(eb0, ea0, ea1, local_grad);
        grad.head(2 * dim) = local_grad.tail(2 * dim);
        grad.segment(2 * dim, dim) = local_grad.head(dim);
        break;

    case EdgeEdgeDistanceType::EA_EB1:
        point_edge_distance_gradient(eb1, ea0, ea1, local_grad);
        grad.head(2 * dim) = local_grad.tail(2 * dim);
        grad.tail(dim) = local_grad.head(dim);
        break;

    case EdgeEdgeDistanceType::EA0_EB:
        point_edge_distance_gradient(ea0, eb0, eb1, local_grad);
        grad.head(dim) = local_grad.head(dim);
        grad.tail(2 * dim) = local_grad.head(2 * dim);
        break;

    case EdgeEdgeDistanceType::EA1_EB:
        point_edge_distance_gradient(ea1, eb0, eb1, local_grad);
        grad.tail(3 * dim) = local_grad;
        break;

    case EdgeEdgeDistanceType::EA_EB:
        line_line_distance_gradient(ea0, ea1, eb0, eb1, grad);
        break;
    }
}

/// @brief Compute the hessian of the distance between a two lines segments.
/// @note The distance is actually squared distance.
/// @param[in] ea0,ea1 The points of the first edge.
/// @param[in] eb0,eb1 The points of the second edge.
/// @param[out] hess The computed hessian.
/// @param[in] project_to_psd True if the hessian should be projected to
///                           positive semi-definite.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1,
    typename DerivedHess>
auto edge_edge_distance_hessian(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1,
    Eigen::MatrixBase<DerivedHess>& hess,
    bool project_to_psd = false)
{
    int dim = ea0.size();
    assert(ea1.size() == dim);
    assert(eb0.size() == dim);
    assert(eb1.size() == dim);

    hess.resize(4 * dim, 4 * dim);
    hess.setZero();

    Eigen::MatrixXd local_hess;
    switch (edge_edge_distance_type(ea0, ea1, eb0, eb1)) {
    case EdgeEdgeDistanceType::EA0_EB0:
        point_point_distance_hessian(ea0, eb0, local_hess, project_to_psd);
        hess.topLeftCorner(dim, dim) = local_hess.topLeftCorner(dim, dim);
        hess.block(0, 2 * dim, dim, dim) = local_hess.topRightCorner(dim, dim);
        hess.block(2 * dim, 0, dim, dim) =
            local_hess.bottomLeftCorner(dim, dim);
        hess.block(2 * dim, 2 * dim, dim, dim) =
            local_hess.bottomRightCorner(dim, dim);
        break;

    case EdgeEdgeDistanceType::EA0_EB1:
        point_point_distance_hessian(ea0, eb1, local_hess, project_to_psd);
        hess.topLeftCorner(dim, dim) = local_hess.topLeftCorner(dim, dim);
        hess.topRightCorner(dim, dim) = local_hess.topRightCorner(dim, dim);
        hess.bottomLeftCorner(dim, dim) = local_hess.bottomLeftCorner(dim, dim);
        hess.bottomRightCorner(dim, dim) =
            local_hess.bottomRightCorner(dim, dim);
        break;

    case EdgeEdgeDistanceType::EA1_EB0:
        point_point_distance_hessian(ea1, eb0, local_hess, project_to_psd);
        hess.block(dim, dim, 2 * dim, 2 * dim) = local_hess;
        break;

    case EdgeEdgeDistanceType::EA1_EB1:
        point_point_distance_hessian(ea1, eb1, local_hess, project_to_psd);
        hess.block(dim, dim, dim, dim) = local_hess.topLeftCorner(dim, dim);
        hess.block(dim, 3 * dim, dim, dim) =
            local_hess.topRightCorner(dim, dim);
        hess.block(3 * dim, dim, dim, dim) =
            local_hess.bottomLeftCorner(dim, dim);
        hess.bottomRightCorner(dim, dim) =
            local_hess.bottomRightCorner(dim, dim);
        break;

    case EdgeEdgeDistanceType::EA_EB0:
        point_edge_distance_hessian(eb0, ea0, ea1, local_hess, project_to_psd);
        hess.topLeftCorner(2 * dim, 2 * dim) =
            local_hess.bottomRightCorner(2 * dim, 2 * dim);
        hess.block(0, 2 * dim, 2 * dim, dim) =
            local_hess.bottomLeftCorner(2 * dim, dim);
        hess.block(2 * dim, 0, dim, 2 * dim) =
            local_hess.topRightCorner(dim, 2 * dim);
        hess.block(2 * dim, 2 * dim, dim, dim) =
            local_hess.topLeftCorner(dim, dim);
        break;

    case EdgeEdgeDistanceType::EA_EB1:
        point_edge_distance_hessian(eb1, ea0, ea1, local_hess, project_to_psd);
        hess.topLeftCorner(2 * dim, 2 * dim) =
            local_hess.bottomRightCorner(2 * dim, 2 * dim);
        hess.topRightCorner(2 * dim, dim) =
            local_hess.bottomLeftCorner(2 * dim, dim);
        hess.bottomLeftCorner(dim, 2 * dim) =
            local_hess.topRightCorner(dim, 2 * dim);
        hess.bottomRightCorner(dim, dim) = local_hess.topLeftCorner(dim, dim);
        break;

    case EdgeEdgeDistanceType::EA0_EB:
        point_edge_distance_hessian(ea0, eb0, eb1, local_hess, project_to_psd);
        hess.topLeftCorner(dim, dim) = local_hess.topLeftCorner(dim, dim);
        hess.topRightCorner(dim, 2 * dim) =
            local_hess.topRightCorner(dim, 2 * dim);
        hess.bottomLeftCorner(2 * dim, dim) =
            local_hess.bottomLeftCorner(2 * dim, dim);
        hess.bottomRightCorner(2 * dim, 2 * dim) =
            local_hess.bottomRightCorner(2 * dim, 2 * dim);
        break;

    case EdgeEdgeDistanceType::EA1_EB:
        point_edge_distance_hessian(ea1, eb0, eb1, local_hess, project_to_psd);
        hess.bottomRightCorner(3 * dim, 3 * dim) = local_hess;
        break;

    case EdgeEdgeDistanceType::EA_EB:
        line_line_distance_hessian(ea0, ea1, eb0, eb1, hess, project_to_psd);
        break;
    }
}

} // namespace ipc
