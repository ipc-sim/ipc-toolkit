#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Compute the distance between a point and edge in 2D or 3D.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0 The first vertex of the edge.
/// @param e1 The second vertex of the edge.
/// @param dtype The point edge distance type to compute.
/// @return The distance between the point and edge.
template <typename T>
T point_edge_distance(
    Eigen::ConstRef<VectorMax3<T>> p,
    Eigen::ConstRef<VectorMax3<T>> e0,
    Eigen::ConstRef<VectorMax3<T>> e1,
    PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);

/// @brief Compute the gradient of the distance between a point and edge.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0 The first vertex of the edge.
/// @param e1 The second vertex of the edge.
/// @param dtype The point edge distance type to compute.
/// @return The gradient of the distance wrt p, e0, and e1.
template <typename T>
VectorMax9<T> point_edge_distance_gradient(
    Eigen::ConstRef<VectorMax3<T>> p,
    Eigen::ConstRef<VectorMax3<T>> e0,
    Eigen::ConstRef<VectorMax3<T>> e1,
    PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);

/// @brief Compute the hessian of the distance between a point and edge.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0 The first vertex of the edge.
/// @param e1 The second vertex of the edge.
/// @param dtype The point edge distance type to compute.
/// @return The hessian of the distance wrt p, e0, and e1.
template <typename T>
MatrixMax9<T> point_edge_distance_hessian(
    Eigen::ConstRef<VectorMax3<T>> p,
    Eigen::ConstRef<VectorMax3<T>> e0,
    Eigen::ConstRef<VectorMax3<T>> e1,
    PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);

// --- EigenExpression wrappers ---

template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_edge_distance(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1,
    PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO) ->
    typename DerivedP::Scalar
{
    using T = typename DerivedP::Scalar;
    return point_edge_distance(
        Eigen::Ref<const VectorMax3<T>>(p), Eigen::Ref<const VectorMax3<T>>(e0),
        Eigen::Ref<const VectorMax3<T>>(e1), dtype);
}

template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_edge_distance_gradient(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1,
    PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO)
    -> VectorMax9<typename DerivedP::Scalar>
{
    using T = typename DerivedP::Scalar;
    return point_edge_distance_gradient(
        Eigen::Ref<const VectorMax3<T>>(p), Eigen::Ref<const VectorMax3<T>>(e0),
        Eigen::Ref<const VectorMax3<T>>(e1), dtype);
}

template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_edge_distance_hessian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1,
    PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO)
    -> MatrixMax9<typename DerivedP::Scalar>
{
    using T = typename DerivedP::Scalar;
    return point_edge_distance_hessian(
        Eigen::Ref<const VectorMax3<T>>(p), Eigen::Ref<const VectorMax3<T>>(e0),
        Eigen::Ref<const VectorMax3<T>>(e1), dtype);
}

} // namespace ipc
