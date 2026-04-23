#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Compute the distance between two lines segments in 3D.
/// @note The distance is actually squared distance.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @param dtype The point edge distance type to compute.
/// @return The distance between the two edges.
template <typename T>
T edge_edge_distance(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1,
    EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO);

/// @brief Compute the gradient of the distance between a two lines segments.
/// @note The distance is actually squared distance.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @param dtype The point edge distance type to compute.
/// @return The gradient of the distance wrt ea0, ea1, eb0, and eb1.
template <typename T>
Eigen::Vector<T, 12> edge_edge_distance_gradient(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1,
    EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO);

/// @brief Compute the hessian of the distance between a two lines segments.
/// @note The distance is actually squared distance.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @param dtype The point edge distance type to compute.
/// @return The hessian of the distance wrt ea0, ea1, eb0, and eb1.
template <typename T>
Eigen::Matrix<T, 12, 12> edge_edge_distance_hessian(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1,
    EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO);

// --- EigenExpression wrappers ---

template <
    EigenExpression DerivedEa0,
    EigenExpression DerivedEa1,
    EigenExpression DerivedEb0,
    EigenExpression DerivedEb1>
inline auto edge_edge_distance(
    const Eigen::MatrixBase<DerivedEa0>& ea0,
    const Eigen::MatrixBase<DerivedEa1>& ea1,
    const Eigen::MatrixBase<DerivedEb0>& eb0,
    const Eigen::MatrixBase<DerivedEb1>& eb1,
    EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO) ->
    typename DerivedEa0::Scalar
{
    using T = typename DerivedEa0::Scalar;
    return edge_edge_distance(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1), dtype);
}

template <
    EigenExpression DerivedEa0,
    EigenExpression DerivedEa1,
    EigenExpression DerivedEb0,
    EigenExpression DerivedEb1>
inline auto edge_edge_distance_gradient(
    const Eigen::MatrixBase<DerivedEa0>& ea0,
    const Eigen::MatrixBase<DerivedEa1>& ea1,
    const Eigen::MatrixBase<DerivedEb0>& eb0,
    const Eigen::MatrixBase<DerivedEb1>& eb1,
    EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO)
    -> Eigen::Vector<typename DerivedEa0::Scalar, 12>
{
    using T = typename DerivedEa0::Scalar;
    return edge_edge_distance_gradient(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1), dtype);
}

template <
    EigenExpression DerivedEa0,
    EigenExpression DerivedEa1,
    EigenExpression DerivedEb0,
    EigenExpression DerivedEb1>
inline auto edge_edge_distance_hessian(
    const Eigen::MatrixBase<DerivedEa0>& ea0,
    const Eigen::MatrixBase<DerivedEa1>& ea1,
    const Eigen::MatrixBase<DerivedEb0>& eb0,
    const Eigen::MatrixBase<DerivedEb1>& eb1,
    EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO)
    -> Eigen::Matrix<typename DerivedEa0::Scalar, 12, 12>
{
    using T = typename DerivedEa0::Scalar;
    return edge_edge_distance_hessian(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1), dtype);
}

} // namespace ipc
