#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Closest pair between a point and point.
enum class PointPointDistanceType : uint8_t {
    P_P = 0, ///< The points are closest to each other.
    AUTO = 0 ///< Automatically determine the closest pair. (Same as P_P)
};

/// @brief Closest pair between a point and edge.
enum class PointEdgeDistanceType : uint8_t {
    P_E0, ///< The point is closest to edge vertex zero.
    P_E1, ///< The point is closest to edge vertex one.
    P_E,  ///< The point is closest to the interior of the edge.
    AUTO  ///< Automatically determine the closest pair.
};

/// @brief Closest pair between a point and triangle.
enum class PointTriangleDistanceType : uint8_t {
    P_T0, ///< The point is closest to triangle vertex zero.
    P_T1, ///< The point is closest to triangle vertex one.
    P_T2, ///< The point is closest to triangle vertex two.
    P_E0, ///< The point is closest to triangle edge zero (vertex zero to one).
    P_E1, ///< The point is closest to triangle edge one (vertex one to two).
    P_E2, ///< The point is closest to triangle edge two (vertex two to zero).
    P_T,  ///< The point is closest to the interior of the triangle.
    AUTO  ///< Automatically determine the closest pair.
};

/// @brief Closest pair between two edges.
enum class EdgeEdgeDistanceType : uint8_t {
    /// The edges are closest at vertex 0 of edge A and 0 of edge B.
    EA0_EB0,
    /// The edges are closest at vertex 0 of edge A and 1 of edge B.
    EA0_EB1,
    /// The edges are closest at vertex 1 of edge A and 0 of edge B.
    EA1_EB0,
    /// The edges are closest at vertex 1 of edge A and 1 of edge B.
    EA1_EB1,
    /// The edges are closest at the interior of edge A and vertex 0 of edge B.
    EA_EB0,
    /// The edges are closest at the interior of edge A and vertex 1 of edge B.
    EA_EB1,
    /// The edges are closest at vertex 0 of edge A and the interior of edge B.
    EA0_EB,
    /// The edges are closest at vertex 1 of edge A and the interior of edge B.
    EA1_EB,
    /// The edges are closest at an interior point of edge A and B.
    EA_EB,
    /// Automatically determine the closest pair.
    AUTO
};

/// @brief Determine the closest pair between a point and edge.
/// @param p The point.
/// @param e0 The first vertex of the edge.
/// @param e1 The second vertex of the edge.
/// @return The distance type of the point-edge pair.
template <typename T>
PointEdgeDistanceType point_edge_distance_type(
    Eigen::ConstRef<VectorMax3<T>> p,
    Eigen::ConstRef<VectorMax3<T>> e0,
    Eigen::ConstRef<VectorMax3<T>> e1);

/// @brief Determine the closest pair between a point and triangle.
/// @param p The point.
/// @param t0 The first vertex of the triangle.
/// @param t1 The second vertex of the triangle.
/// @param t2 The third vertex of the triangle.
/// @return The distance type of the point-triangle pair.
template <typename T>
PointTriangleDistanceType point_triangle_distance_type(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> t0,
    Eigen::ConstRef<Eigen::Vector3<T>> t1,
    Eigen::ConstRef<Eigen::Vector3<T>> t2);

/// @brief Determine the closest pair between two edges.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @return The distance type of the edge-edge pair.
template <typename T>
EdgeEdgeDistanceType edge_edge_distance_type(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1);

/// @brief Determine the closest pair between two parallel edges.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @return The distance type of the edge-edge pair.
template <typename T>
EdgeEdgeDistanceType edge_edge_parallel_distance_type(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1);

// --- EigenExpression wrappers ---

template <EigenExpression D0, EigenExpression D1, EigenExpression D2>
inline PointEdgeDistanceType point_edge_distance_type(
    const Eigen::MatrixBase<D0>& p,
    const Eigen::MatrixBase<D1>& e0,
    const Eigen::MatrixBase<D2>& e1)
{
    using T = typename D0::Scalar;
    return point_edge_distance_type(
        Eigen::Ref<const VectorMax3<T>>(p), Eigen::Ref<const VectorMax3<T>>(e0),
        Eigen::Ref<const VectorMax3<T>>(e1));
}

template <
    EigenExpression D0,
    EigenExpression D1,
    EigenExpression D2,
    EigenExpression D3>
inline PointTriangleDistanceType point_triangle_distance_type(
    const Eigen::MatrixBase<D0>& p,
    const Eigen::MatrixBase<D1>& t0,
    const Eigen::MatrixBase<D2>& t1,
    const Eigen::MatrixBase<D3>& t2)
{
    using T = typename D0::Scalar;
    return point_triangle_distance_type(
        Eigen::Ref<const Eigen::Vector3<T>>(p),
        Eigen::Ref<const Eigen::Vector3<T>>(t0),
        Eigen::Ref<const Eigen::Vector3<T>>(t1),
        Eigen::Ref<const Eigen::Vector3<T>>(t2));
}

template <
    EigenExpression D0,
    EigenExpression D1,
    EigenExpression D2,
    EigenExpression D3>
inline EdgeEdgeDistanceType edge_edge_distance_type(
    const Eigen::MatrixBase<D0>& ea0,
    const Eigen::MatrixBase<D1>& ea1,
    const Eigen::MatrixBase<D2>& eb0,
    const Eigen::MatrixBase<D3>& eb1)
{
    using T = typename D0::Scalar;
    return edge_edge_distance_type(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1));
}

template <
    EigenExpression D0,
    EigenExpression D1,
    EigenExpression D2,
    EigenExpression D3>
inline EdgeEdgeDistanceType edge_edge_parallel_distance_type(
    const Eigen::MatrixBase<D0>& ea0,
    const Eigen::MatrixBase<D1>& ea1,
    const Eigen::MatrixBase<D2>& eb0,
    const Eigen::MatrixBase<D3>& eb1)
{
    using T = typename D0::Scalar;
    return edge_edge_parallel_distance_type(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1));
}

} // namespace ipc
