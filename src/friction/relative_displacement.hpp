#pragma once

#include <Eigen/Core>

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

// Point - Triangle

/// Compute the relative displacement of the point to the triangle.
template <
    typename DerivedDP,
    typename DerivedDT0,
    typename DerivedDT1,
    typename DerivedDT2,
    typename DerivedBarycentricCoordinates>
inline auto point_triangle_relative_displacement(
    const Eigen::MatrixBase<DerivedDP>& dp,
    const Eigen::MatrixBase<DerivedDT0>& dt0,
    const Eigen::MatrixBase<DerivedDT1>& dt1,
    const Eigen::MatrixBase<DerivedDT2>& dt2,
    const Eigen::MatrixBase<DerivedBarycentricCoordinates>&
        barycentric_coordinates)
{
    // Compute the displacement of the closest point and subtract it from the
    // points displacement.
    return dp
        - (dt0 + barycentric_coordinates[0] * (dt1 - dt0)
           + barycentric_coordinates[1] * (dt2 - dt0));
}

template <
    typename DerivedDisp,
    typename DerivedBasis,
    typename DerivedBeta,
    typename T = typename DerivedDisp::Scalar>
inline Eigen::Matrix<T, 4, 3> point_triangle_relative_mesh_displacements(
    const Eigen::MatrixBase<DerivedDisp>& tangent_relative_displacement,
    const Eigen::MatrixBase<DerivedBasis>& basis,
    const Eigen::MatrixBase<DerivedBeta>& beta)
{
    Eigen::Matrix<T, 4, 3> mesh_displacements;
    mesh_displacements.row(0) = basis * tangent_relative_displacement;
    mesh_displacements.row(1) =
        (-1 + beta[0] + beta[1]) * mesh_displacements.row(0);
    mesh_displacements.row(2) = -beta[0] * mesh_displacements.row(0);
    mesh_displacements.row(3) = -beta[1] * mesh_displacements.row(0);
    return mesh_displacements;
}

template <typename DerivedBasis, typename DerivedBeta, typename DerivedTT>
inline void point_triangle_TT(
    const Eigen::MatrixBase<DerivedBasis>& basis,
    const Eigen::MatrixBase<DerivedBeta>& beta,
    Eigen::MatrixBase<DerivedTT>& TT)
{
    TT.template block<2, 3>(0, 0) = basis.transpose();
    TT.template block<2, 3>(0, 3) =
        (-1 + beta[0] + beta[1]) * basis.transpose();
    TT.template block<2, 3>(0, 6) = -beta[0] * basis.transpose();
    TT.template block<2, 3>(0, 9) = -beta[1] * basis.transpose();
}

// Edge - Edge

/// Compute the relative displacement of the edges.
template <
    typename DerivedDEA0,
    typename DerivedDEA1,
    typename DerivedDEB0,
    typename DerivedDEB1,
    typename DerivedBarycentricCoordinates>
inline auto edge_edge_relative_displacement(
    const Eigen::MatrixBase<DerivedDEA0>& dea0,
    const Eigen::MatrixBase<DerivedDEA1>& dea1,
    const Eigen::MatrixBase<DerivedDEB0>& deb0,
    const Eigen::MatrixBase<DerivedDEB1>& deb1,
    const Eigen::MatrixBase<DerivedBarycentricCoordinates>&
        barrycentric_coordinates)
{
    // closest_point_a_displacement - closest_point_b_displacement
    return ((dea1 - dea0) * barrycentric_coordinates[0] + dea0)
        - ((deb1 - deb0) * barrycentric_coordinates[1] + deb0);
}

template <
    typename DerivedDisp,
    typename DerivedBasis,
    typename DerivedBeta,
    typename T = typename DerivedDisp::Scalar>
inline Eigen::Matrix<T, 4, 3> edge_edge_relative_mesh_displacements(
    const Eigen::MatrixBase<DerivedDisp>& tangent_relative_displacement,
    const Eigen::MatrixBase<DerivedBasis>& basis,
    const Eigen::MatrixBase<DerivedBeta>& gamma)
{
    Eigen::Vector3<T> rel_disp = basis * tangent_relative_displacement;
    Eigen::Matrix<T, 4, 3> mesh_displacements;
    mesh_displacements.row(0) = (1.0 - gamma[0]) * rel_disp; // dea0
    mesh_displacements.row(1) = gamma[0] * rel_disp;         // dea1
    mesh_displacements.row(2) = (gamma[1] - 1.0) * rel_disp; // deb0
    mesh_displacements.row(3) = -gamma[1] * rel_disp;        // deb1
    return mesh_displacements;
}

template <typename DerivedBasis, typename DerivedGamma, typename DerivedTT>
inline void edge_edge_TT(
    const Eigen::MatrixBase<DerivedBasis>& basis,
    const Eigen::MatrixBase<DerivedGamma>& gamma,
    Eigen::MatrixBase<DerivedTT>& TT)
{
    TT.template block<2, 3>(0, 0) = (1.0 - gamma[0]) * basis.transpose();
    TT.template block<2, 3>(0, 3) = gamma[0] * basis.transpose();
    TT.template block<2, 3>(0, 6) = (gamma[1] - 1.0) * basis.transpose();
    TT.template block<2, 3>(0, 9) = -gamma[1] * basis.transpose();
}

// Point - Edge

template <
    typename DerivedDP,
    typename DerivedDE0,
    typename DerivedDE1,
    typename T>
inline auto point_edge_relative_displacement(
    const Eigen::MatrixBase<DerivedDP>& dp,
    const Eigen::MatrixBase<DerivedDE0>& de0,
    const Eigen::MatrixBase<DerivedDE1>& de1,
    const T& alpha)
{
    return dp - ((de1 - de0) * alpha + de0);
}

template <typename DerivedDisp, typename DerivedBasis, typename T>
inline Eigen::Matrix3<T> point_edge_relative_mesh_displacement(
    const Eigen::MatrixBase<DerivedDisp>& tangent_relative_displacement,
    const Eigen::MatrixBase<DerivedBasis>& basis,
    const T& alpha)
{
    Eigen::Matrix3<T> mesh_displacements;
    mesh_displacements.row(0) = basis * tangent_relative_displacement;
    mesh_displacements.row(1) = (alpha - 1.0) * mesh_displacements.row(0);
    mesh_displacements.row(2) = -alpha * mesh_displacements.row(0);
    return mesh_displacements;
}

template <typename DerivedBasis, typename T, typename DerivedTT>
inline void point_edge_TT(
    const Eigen::MatrixBase<DerivedBasis>& basis,
    const T& alpha,
    Eigen::MatrixBase<DerivedTT>& TT)
{
    TT.template block<2, 3>(0, 0) = basis.transpose();
    TT.template block<2, 3>(0, 3) = (alpha - 1.0) * basis.transpose();
    TT.template block<2, 3>(0, 6) = -alpha * basis.transpose();
}

// Point - Point

template <typename DerivedDP0, typename DerivedDP1>
inline auto point_point_relative_displacement(
    const Eigen::MatrixBase<DerivedDP0>& dp0,
    const Eigen::MatrixBase<DerivedDP1>& dp1)
{
    return dp0 - dp1;
}

template <typename T>
inline void point_point_relative_mesh_displacement(
    const Eigen::Vector2<T>& tangent_relative_displacement,
    const Eigen::Matrix<T, 3, 2>& basis,
    Eigen::Matrix<T, 2, 3>& mesh_displacements)
{
    mesh_displacements.row(0) = basis * tangent_relative_displacement;
    mesh_displacements.row(1) = -mesh_displacements.row(0);
}

template <typename DerivedBasis, typename DerivedTT>
inline void point_point_TT(
    const Eigen::MatrixBase<DerivedBasis>& basis,
    Eigen::MatrixBase<DerivedTT>& TT)
{
    TT.template block<2, 3>(0, 0) = basis.transpose();
    TT.template block<2, 3>(0, 3) = -basis.transpose();
}

} // namespace ipc
