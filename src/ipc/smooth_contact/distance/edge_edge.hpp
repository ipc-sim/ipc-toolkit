#pragma once

#include "point_edge.hpp"

namespace ipc {

template <typename T>
T line_line_sqr_distance(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1);

template <typename T>
T edge_edge_sqr_distance(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1,
    EdgeEdgeDistanceType dtype);

template <typename T>
Eigen::Vector3<T> line_line_closest_point_direction(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1);

std::tuple<Eigen::Vector3d, Eigen::Matrix<double, 3, 12>>
line_line_closest_point_direction_gradient(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1);

std::tuple<
    Eigen::Vector3d,
    Eigen::Matrix<double, 3, 12>,
    std::array<Matrix12d, 3>>
line_line_closest_point_direction_hessian(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1);

template <typename T>
Eigen::Matrix<T, 3, 2> line_line_closest_point_pairs(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1);

std::tuple<Vector6d, Eigen::Matrix<double, 6, 12>>
line_line_closest_point_pairs_gradient(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1);

std::tuple<Vector6d, Eigen::Matrix<double, 6, 12>, std::array<Matrix12d, 6>>
line_line_closest_point_pairs_hessian(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1);

/// @brief Computes the direction of the closest point pair
/// @param ea0 Vertex 0 of edge 0
/// @param ea1 Vertex 1 of edge 0
/// @param eb0 Vertex 0 of edge 1
/// @param eb1 Vertex 1 of edge 1
/// @param dtype Edge-edge distance type
/// @return Difference of the pair of closest point, pointing from edge 0 to edge 1
template <typename T>
Eigen::Vector3<T> edge_edge_closest_point_direction(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1,
    EdgeEdgeDistanceType dtype);

/// @brief Computes the position of two closest points on two edges
/// @param ea0 Vertex 0 of edge 0
/// @param ea1 Vertex 1 of edge 0
/// @param eb0 Vertex 0 of edge 1
/// @param eb1 Vertex 1 of edge 1
/// @param dtype Edge-edge distance type
template <typename T>
Eigen::Matrix<T, 3, 2> edge_edge_closest_point_pairs(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1,
    EdgeEdgeDistanceType dtype);
} // namespace ipc

#include "edge_edge.tpp"
