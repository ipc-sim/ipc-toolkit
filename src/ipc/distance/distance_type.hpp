#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Closest pair between a point and edge.
enum class PointEdgeDistanceType {
    P_E0, ///< The point is closest to edge vertex zero.
    P_E1, ///< The point is closest to edge vertex one.
    P_E,  ///< The point is closest to the interior of the edge.
    AUTO  ///< Automatically determine the closest pair.
};

/// @brief Closest pair between a point and triangle.
enum class PointTriangleDistanceType {
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
enum class EdgeEdgeDistanceType {
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
PointEdgeDistanceType point_edge_distance_type(
    const Eigen::Ref<const VectorMax3d>& p,
    const Eigen::Ref<const VectorMax3d>& e0,
    const Eigen::Ref<const VectorMax3d>& e1);

/// @brief Determine the closest pair between a point and triangle.
/// @param p The point.
/// @param t0 The first vertex of the triangle.
/// @param t1 The second vertex of the triangle.
/// @param t2 The third vertex of the triangle.
/// @return The distance type of the point-triangle pair.
PointTriangleDistanceType point_triangle_distance_type(
    const Eigen::Ref<const Eigen::Vector3d>& p,
    const Eigen::Ref<const Eigen::Vector3d>& t0,
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2);

/// @brief Determine the closest pair between two edges.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @return The distance type of the edge-edge pair.
EdgeEdgeDistanceType edge_edge_distance_type(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1);

/// @brief Determine the closest pair between two parallel edges.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @return The distance type of the edge-edge pair.
EdgeEdgeDistanceType edge_edge_parallel_distance_type(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1);

} // namespace ipc
