#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

namespace detail {
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
        EdgeEdgeDistanceType dtype);

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
        EdgeEdgeDistanceType dtype);

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
        EdgeEdgeDistanceType dtype);
} // namespace detail

// --- EigenExpression wrappers ---

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto edge_edge_distance(
    const DerivedEA0& ea0,
    const DerivedEA1& ea1,
    const DerivedEB0& eb0,
    const DerivedEB1& eb1,
    EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedEA0, DerivedEA1, DerivedEB0, DerivedEB1);
    using T = typename DerivedEA0::Scalar;
    return detail::edge_edge_distance<T>(ea0, ea1, eb0, eb1, dtype);
}

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto edge_edge_distance_gradient(
    const DerivedEA0& ea0,
    const DerivedEA1& ea1,
    const DerivedEB0& eb0,
    const DerivedEB1& eb1,
    EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedEA0, DerivedEA1, DerivedEB0, DerivedEB1);
    using T = typename DerivedEA0::Scalar;
    return detail::edge_edge_distance_gradient<T>(ea0, ea1, eb0, eb1, dtype);
}

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto edge_edge_distance_hessian(
    const DerivedEA0& ea0,
    const DerivedEA1& ea1,
    const DerivedEB0& eb0,
    const DerivedEB1& eb1,
    EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedEA0, DerivedEA1, DerivedEB0, DerivedEB1);
    using T = typename DerivedEA0::Scalar;
    return detail::edge_edge_distance_hessian<T>(ea0, ea1, eb0, eb1, dtype);
}

} // namespace ipc
