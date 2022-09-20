#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/line_line.hpp>

#include <stdexcept> // std::invalid_argument

namespace ipc {

/// @brief Compute the distance between a two lines segments in 3D.
/// @note The distance is actually squared distance.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @return The distance between the two edges.
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
    assert(ea0.size() == 3);
    assert(ea1.size() == 3);
    assert(eb0.size() == 3);
    assert(eb1.size() == 3);

    return edge_edge_distance(
        ea0, ea1, eb0, eb1, edge_edge_distance_type(ea0, ea1, eb0, eb1));
}

/// @brief Compute the distance between a two lines segments in 3D.
/// @note The distance is actually squared distance.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @param dtype The point edge distance type to compute.
/// @return The distance between the two edges.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
auto edge_edge_distance(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1,
    const EdgeEdgeDistanceType dtype)
{
    assert(ea0.size() == 3);
    assert(ea1.size() == 3);
    assert(eb0.size() == 3);
    assert(eb1.size() == 3);

    switch (dtype) {
    case EdgeEdgeDistanceType::EA0_EB0:
        return point_point_distance(ea0, eb0);

    case EdgeEdgeDistanceType::EA0_EB1:
        return point_point_distance(ea0, eb1);

    case EdgeEdgeDistanceType::EA1_EB0:
        return point_point_distance(ea1, eb0);

    case EdgeEdgeDistanceType::EA1_EB1:
        return point_point_distance(ea1, eb1);

    case EdgeEdgeDistanceType::EA_EB0:
        return point_line_distance(eb0, ea0, ea1);

    case EdgeEdgeDistanceType::EA_EB1:
        return point_line_distance(eb1, ea0, ea1);

    case EdgeEdgeDistanceType::EA0_EB:
        return point_line_distance(ea0, eb0, eb1);

    case EdgeEdgeDistanceType::EA1_EB:
        return point_line_distance(ea1, eb0, eb1);

    case EdgeEdgeDistanceType::EA_EB:
        return line_line_distance(ea0, ea1, eb0, eb1);

    default:
        throw std::invalid_argument(
            "Invalid distance type for edge-edge distance!");
    }
}

/// @brief Compute the gradient of the distance between a two lines segments.
/// @note The distance is actually squared distance.
/// @param[in] ea0 The first vertex of the first edge.
/// @param[in] ea1 The second vertex of the first edge.
/// @param[in] eb0 The first vertex of the second edge.
/// @param[in] eb1 The second vertex of the second edge.
/// @param[out] grad The gradient of the distance wrt ea0, ea1, eb0, and eb1.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1,
    typename DerivedGrad>
void edge_edge_distance_gradient(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1,
    Eigen::PlainObjectBase<DerivedGrad>& grad)
{
    assert(ea0.size() == 3);
    assert(ea1.size() == 3);
    assert(eb0.size() == 3);
    assert(eb1.size() == 3);

    return edge_edge_distance_gradient(
        ea0, ea1, eb0, eb1, edge_edge_distance_type(ea0, ea1, eb0, eb1), grad);
}

/// @brief Compute the gradient of the distance between a two lines segments.
/// @note The distance is actually squared distance.
/// @param[in] ea0 The first vertex of the first edge.
/// @param[in] ea1 The second vertex of the first edge.
/// @param[in] eb0 The first vertex of the second edge.
/// @param[in] eb1 The second vertex of the second edge.
/// @param[in] dtype The point edge distance type to compute.
/// @param[out] grad The gradient of the distance wrt ea0, ea1, eb0, and eb1.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1,
    typename DerivedGrad>
void edge_edge_distance_gradient(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1,
    const EdgeEdgeDistanceType dtype,
    Eigen::PlainObjectBase<DerivedGrad>& grad)
{
    int dim = ea0.size();
    assert(ea1.size() == dim);
    assert(eb0.size() == dim);
    assert(eb1.size() == dim);

    grad.resize(4 * dim);
    grad.setZero();

    VectorMax9<typename DerivedGrad::Scalar> local_grad;
    switch (dtype) {
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
        point_line_distance_gradient(eb0, ea0, ea1, local_grad);
        grad.head(2 * dim) = local_grad.tail(2 * dim);
        grad.segment(2 * dim, dim) = local_grad.head(dim);
        break;

    case EdgeEdgeDistanceType::EA_EB1:
        point_line_distance_gradient(eb1, ea0, ea1, local_grad);
        grad.head(2 * dim) = local_grad.tail(2 * dim);
        grad.tail(dim) = local_grad.head(dim);
        break;

    case EdgeEdgeDistanceType::EA0_EB:
        point_line_distance_gradient(ea0, eb0, eb1, local_grad);
        grad.head(dim) = local_grad.head(dim);
        grad.tail(2 * dim) = local_grad.tail(2 * dim);
        break;

    case EdgeEdgeDistanceType::EA1_EB:
        point_line_distance_gradient(ea1, eb0, eb1, local_grad);
        grad.tail(3 * dim) = local_grad;
        break;

    case EdgeEdgeDistanceType::EA_EB:
        line_line_distance_gradient(ea0, ea1, eb0, eb1, grad);
        break;

    default:
        throw std::invalid_argument(
            "Invalid distance type for edge-edge distance gradient!");
    }
}

/// @brief Compute the hessian of the distance between a two lines segments.
/// @note The distance is actually squared distance.
/// @param[in] ea0 The first vertex of the first edge.
/// @param[in] ea1 The second vertex of the first edge.
/// @param[in] eb0 The first vertex of the second edge.
/// @param[in] eb1 The second vertex of the second edge.
/// @param[out] hess The hessian of the distance wrt ea0, ea1, eb0, and eb1.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1,
    typename DerivedHess>
void edge_edge_distance_hessian(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1,
    Eigen::PlainObjectBase<DerivedHess>& hess)
{
    assert(ea0.size() == 3);
    assert(ea1.size() == 3);
    assert(eb0.size() == 3);
    assert(eb1.size() == 3);

    return edge_edge_distance_hessian(
        ea0, ea1, eb0, eb1, edge_edge_distance_type(ea0, ea1, eb0, eb1), hess);
}

/// @brief Compute the hessian of the distance between a two lines segments.
/// @note The distance is actually squared distance.
/// @param[in] ea0 The first vertex of the first edge.
/// @param[in] ea1 The second vertex of the first edge.
/// @param[in] eb0 The first vertex of the second edge.
/// @param[in] eb1 The second vertex of the second edge.
/// @param[in] dtype The point edge distance type to compute.
/// @param[out] hess The hessian of the distance wrt ea0, ea1, eb0, and eb1.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1,
    typename DerivedHess>
void edge_edge_distance_hessian(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1,
    const EdgeEdgeDistanceType dtype,
    Eigen::PlainObjectBase<DerivedHess>& hess)
{
    int dim = ea0.size();
    assert(ea1.size() == dim);
    assert(eb0.size() == dim);
    assert(eb1.size() == dim);

    hess.resize(4 * dim, 4 * dim);
    hess.setZero();

    MatrixMax9<typename DerivedHess::Scalar> local_hess;
    switch (dtype) {
    case EdgeEdgeDistanceType::EA0_EB0:
        point_point_distance_hessian(ea0, eb0, local_hess);
        hess.topLeftCorner(dim, dim) = local_hess.topLeftCorner(dim, dim);
        hess.block(0, 2 * dim, dim, dim) = local_hess.topRightCorner(dim, dim);
        hess.block(2 * dim, 0, dim, dim) =
            local_hess.bottomLeftCorner(dim, dim);
        hess.block(2 * dim, 2 * dim, dim, dim) =
            local_hess.bottomRightCorner(dim, dim);
        break;

    case EdgeEdgeDistanceType::EA0_EB1:
        point_point_distance_hessian(ea0, eb1, local_hess);
        hess.topLeftCorner(dim, dim) = local_hess.topLeftCorner(dim, dim);
        hess.topRightCorner(dim, dim) = local_hess.topRightCorner(dim, dim);
        hess.bottomLeftCorner(dim, dim) = local_hess.bottomLeftCorner(dim, dim);
        hess.bottomRightCorner(dim, dim) =
            local_hess.bottomRightCorner(dim, dim);
        break;

    case EdgeEdgeDistanceType::EA1_EB0:
        point_point_distance_hessian(ea1, eb0, local_hess);
        hess.block(dim, dim, 2 * dim, 2 * dim) = local_hess;
        break;

    case EdgeEdgeDistanceType::EA1_EB1:
        point_point_distance_hessian(ea1, eb1, local_hess);
        hess.block(dim, dim, dim, dim) = local_hess.topLeftCorner(dim, dim);
        hess.block(dim, 3 * dim, dim, dim) =
            local_hess.topRightCorner(dim, dim);
        hess.block(3 * dim, dim, dim, dim) =
            local_hess.bottomLeftCorner(dim, dim);
        hess.bottomRightCorner(dim, dim) =
            local_hess.bottomRightCorner(dim, dim);
        break;

    case EdgeEdgeDistanceType::EA_EB0:
        point_line_distance_hessian(eb0, ea0, ea1, local_hess);
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
        point_line_distance_hessian(eb1, ea0, ea1, local_hess);
        hess.topLeftCorner(2 * dim, 2 * dim) =
            local_hess.bottomRightCorner(2 * dim, 2 * dim);
        hess.topRightCorner(2 * dim, dim) =
            local_hess.bottomLeftCorner(2 * dim, dim);
        hess.bottomLeftCorner(dim, 2 * dim) =
            local_hess.topRightCorner(dim, 2 * dim);
        hess.bottomRightCorner(dim, dim) = local_hess.topLeftCorner(dim, dim);
        break;

    case EdgeEdgeDistanceType::EA0_EB:
        point_line_distance_hessian(ea0, eb0, eb1, local_hess);
        hess.topLeftCorner(dim, dim) = local_hess.topLeftCorner(dim, dim);
        hess.topRightCorner(dim, 2 * dim) =
            local_hess.topRightCorner(dim, 2 * dim);
        hess.bottomLeftCorner(2 * dim, dim) =
            local_hess.bottomLeftCorner(2 * dim, dim);
        hess.bottomRightCorner(2 * dim, 2 * dim) =
            local_hess.bottomRightCorner(2 * dim, 2 * dim);
        break;

    case EdgeEdgeDistanceType::EA1_EB:
        point_line_distance_hessian(ea1, eb0, eb1, local_hess);
        hess.bottomRightCorner(3 * dim, 3 * dim) = local_hess;
        break;

    case EdgeEdgeDistanceType::EA_EB:
        line_line_distance_hessian(ea0, ea1, eb0, eb1, hess);
        break;

    default:
        throw std::invalid_argument(
            "Invalid distance type for edge-edge distance hessian!");
    }
}

} // namespace ipc
