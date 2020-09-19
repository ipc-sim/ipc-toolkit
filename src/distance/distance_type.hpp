#pragma once

#include <Eigen/Core>

namespace ipc {

/// @brief Closest pair between a point and edge.
enum class PointEdgeDistanceType { P_E0, P_E1, P_E };

/// @brief Closest pair between a point and triangle.
enum class PointTriangleDistanceType {
    P_T0,
    P_T1,
    P_T2,
    P_E0,
    P_E1,
    P_E2,
    P_T
};

/// @brief Closest pair between two edges.
enum class EdgeEdgeDistanceType {
    EA0_EB0,
    EA0_EB1,
    EA1_EB0,
    EA1_EB1,
    EA_EB0,
    EA_EB1,
    EA0_EB,
    EA1_EB,
    EA_EB
};

/// @brief Determine the closest pair between a point and edge.
/// @param p The point.
/// @param e0,e1 The points of the edge.
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
PointEdgeDistanceType point_edge_distance_type(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1);

/// @brief Determine the closest pair between a point and triangle.
/// @param p The point.
/// @param t0,t1,t2 The points of the triangle.
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
/// @param ea0,ea1 The points of the first edge.
/// @param eb0,eb1 The points of the second edge.
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
