#pragma once

#include <Eigen/Core>

namespace ipc {

/// @brief Closest pair between a point and edge.
enum class PointEdgeDistanceType {
    P_E0, ///< The point is closest to edge vertex zero.
    P_E1, ///< The point is closest to edge vertex one.
    P_E   ///< The point is closest to the interior of the edge.
};

/// @brief Closest pair between a point and triangle.
enum class PointTriangleDistanceType {
    P_T0, ///< The point is closest to triangle vertex zero.
    P_T1, ///< The point is closest to triangle vertex one.
    P_T2, ///< The point is closest to triangle vertex two.
    P_E0, ///< The point is closest to triangle edge zero (vertex zero to one).
    P_E1, ///< The point is closest to triangle edge one (vertex one to two).
    P_E2, ///< The point is closest to triangle edge two (vertex two to zero).
    P_T   ///< The point is closest to the interior of the triangle.
};

/// @brief Closest pair between two edges.
enum class EdgeEdgeDistanceType {
    EA0_EB0, ///< The edges are closest at vertex 0 of edge A and 0 of edge B.
    EA0_EB1, ///< The edges are closest at vertex 0 of edge A and 1 of edge B.
    EA1_EB0, ///< The edges are closest at vertex 1 of edge A and 0 of edge B.
    EA1_EB1, ///< The edges are closest at vertex 1 of edge A and 1 of edge B.
    /// The edges are closest at the interior of edge A and vertex 0 of edge B.
    EA_EB0,
    /// The edges are closest at the interior of edge A and vertex 1 of edge B.
    EA_EB1,
    /// The edges are closest at vertex 0 of edge A and the interior of edge B.
    EA0_EB,
    /// The edges are closest at vertex 1 of edge A and the interior of edge B.
    EA1_EB,
    EA_EB ///< The edges are closest at an interior point of edge A and B.
};

/// @brief Determine the closest pair between a point and edge.
/// @param p The point.
/// @param e0 The first vertex of the edge.
/// @param e1 The second vertex of the edge.
/// @return The distance type of the point-edge pair.
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
PointEdgeDistanceType point_edge_distance_type(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1);

/// @brief Determine the closest pair between a point and triangle.
/// @param p The point.
/// @param t0 The first vertex of the triangle.
/// @param t1 The second vertex of the triangle.
/// @param t2 The third vertex of the triangle.
/// @return The distance type of the point-triangle pair.
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
PointTriangleDistanceType point_triangle_distance_type(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2);

/// @brief Determine the closest pair between two edges.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @return The distance type of the edge-edge pair.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
EdgeEdgeDistanceType edge_edge_distance_type(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1);

} // namespace ipc

#include <ipc/distance/distance_type.tpp>
