#pragma once

#include <Eigen/Core>

namespace ipc {

enum class PointEdgeDistanceType { P_E0, P_E1, P_E };

enum class PointTriangleDistanceType {
    P_T0,
    P_T1,
    P_T2,
    P_E0,
    P_E1,
    P_E2,
    P_T
};

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

/// Determine the closest pair between a point and edge.
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
PointEdgeDistanceType point_edge_distance_type(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1);

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

#include "distance_type.tpp"
