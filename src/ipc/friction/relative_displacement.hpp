#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

///////////////////////////////////////////////////////////////////////////////
// Point - Point

template <typename DerivedDP0, typename DerivedDP1>
inline auto point_point_relative_displacement(
    const Eigen::MatrixBase<DerivedDP0>& dp0,
    const Eigen::MatrixBase<DerivedDP1>& dp1)
{
    return dp0 - dp1;
}

template <typename T = double>
inline MatrixMax<T, 3, 6>
point_point_relative_displacement_matrix(const int dim)
{
    MatrixMax<T, 3, 6> J(dim, 2 * dim);
    J.leftCols(dim) = MatrixMax<T, 3, 3>::Identity(dim, dim);
    J.rightCols(dim) = -MatrixMax<T, 3, 3>::Identity(dim, dim);
    return J;
}

template <typename T = double>
inline MatrixMax<T, 3, 6>
point_point_relative_displacement_matrix_jacobian(const int dim)
{
    return MatrixMax<T, 3, 6>::Zero(dim, 2 * dim);
}

///////////////////////////////////////////////////////////////////////////////
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

template <typename T>
inline MatrixMax<T, 3, 9>
point_edge_relative_displacement_matrix(const int dim, const T& alpha)
{
    MatrixMax<T, 3, 9> J = MatrixMax<T, 3, 9>::Zero(dim, 3 * dim);
    J.leftCols(dim).diagonal().setOnes();
    J.middleCols(dim, dim).diagonal().setConstant(alpha - 1);
    J.rightCols(dim).diagonal().setConstant(-alpha);
    return J;
}

template <typename T>
inline MatrixMax<T, 3, 9>
point_edge_relative_displacement_matrix_jacobian(const int dim, const T& alpha)
{
    MatrixMax<T, 3, 9> J = MatrixMax<T, 3, 9>::Zero(dim, 3 * dim);
    J.middleCols(dim, dim).diagonal().setConstant(1);
    J.rightCols(dim).diagonal().setConstant(-1);
    return J;
}

///////////////////////////////////////////////////////////////////////////////
// Edge - Edge

/// Compute the relative displacement of the edges.
template <
    typename DerivedDEA0,
    typename DerivedDEA1,
    typename DerivedDEB0,
    typename DerivedDEB1,
    typename DerivedCoords>
inline auto edge_edge_relative_displacement(
    const Eigen::MatrixBase<DerivedDEA0>& dea0,
    const Eigen::MatrixBase<DerivedDEA1>& dea1,
    const Eigen::MatrixBase<DerivedDEB0>& deb0,
    const Eigen::MatrixBase<DerivedDEB1>& deb1,
    const Eigen::MatrixBase<DerivedCoords>& coords)
{
    // closest_point_a_displacement - closest_point_b_displacement
    return ((dea1 - dea0) * coords[0] + dea0)
        - ((deb1 - deb0) * coords[1] + deb0);
}

template <typename DerivedCoords, typename T = typename DerivedCoords::Scalar>
inline MatrixMax<T, 3, 12> edge_edge_relative_displacement_matrix(
    const int dim, const Eigen::MatrixBase<DerivedCoords>& coords)
{
    MatrixMax<T, 3, 12> J = MatrixMax<T, 3, 12>::Zero(dim, 4 * dim);
    J.leftCols(dim).diagonal().setConstant(1 - coords[0]);
    J.middleCols(dim, dim).diagonal().setConstant(coords[0]);
    J.middleCols(2 * dim, dim).diagonal().setConstant(coords[1] - 1);
    J.rightCols(dim).diagonal().setConstant(-coords[1]);
    return J;
}

template <typename DerivedCoords, typename T = typename DerivedCoords::Scalar>
inline MatrixMax<T, 6, 12> edge_edge_relative_displacement_matrix_jacobian(
    const int dim, const Eigen::MatrixBase<DerivedCoords>& coords)
{
    MatrixMax<T, 6, 12> J = MatrixMax<T, 6, 12>::Zero(2 * dim, 4 * dim);
    // wrt β₁
    J.block(0, 0, dim, dim).diagonal().setConstant(-1);
    J.block(0, dim, dim, dim).diagonal().setConstant(1);
    // wrt β₂
    J.block(dim, 2 * dim, dim, dim).diagonal().setConstant(1);
    J.block(dim, 3 * dim, dim, dim).diagonal().setConstant(-1);
    return J;
}

///////////////////////////////////////////////////////////////////////////////
// Point - Triangle

/// Compute the relative displacement of the point to the triangle.
template <
    typename DerivedDP,
    typename DerivedDT0,
    typename DerivedDT1,
    typename DerivedDT2,
    typename DerivedCoords>
inline auto point_triangle_relative_displacement(
    const Eigen::MatrixBase<DerivedDP>& dp,
    const Eigen::MatrixBase<DerivedDT0>& dt0,
    const Eigen::MatrixBase<DerivedDT1>& dt1,
    const Eigen::MatrixBase<DerivedDT2>& dt2,
    const Eigen::MatrixBase<DerivedCoords>& coords)
{
    // Compute the displacement of the closest point and subtract it from the
    // points displacement.
    return dp - (dt0 + coords[0] * (dt1 - dt0) + coords[1] * (dt2 - dt0));
}

template <typename DerivedCoords, typename T = typename DerivedCoords::Scalar>
inline MatrixMax<T, 3, 12> point_triangle_relative_displacement_matrix(
    const int dim, const Eigen::MatrixBase<DerivedCoords>& coords)
{
    MatrixMax<T, 3, 12> J = MatrixMax<T, 3, 12>::Zero(dim, 4 * dim);
    J.leftCols(dim).diagonal().setOnes();
    J.middleCols(dim, dim).diagonal().setConstant(coords[0] + coords[1] - 1);
    J.middleCols(2 * dim, dim).diagonal().setConstant(-coords[0]);
    J.rightCols(dim).diagonal().setConstant(-coords[1]);
    return J;
}

template <typename DerivedCoords, typename T = typename DerivedCoords::Scalar>
inline MatrixMax<T, 6, 12> point_triangle_relative_displacement_matrix_jacobian(
    const int dim, const Eigen::MatrixBase<DerivedCoords>& coords)
{
    MatrixMax<T, 6, 12> J = MatrixMax<T, 6, 12>::Zero(2 * dim, 4 * dim);
    // wrt β₁
    J.block(0, dim, dim, dim).diagonal().setConstant(1);
    J.block(0, 2 * dim, dim, dim).diagonal().setConstant(-1);
    // wrt β₂
    J.block(dim, dim, dim, dim).diagonal().setConstant(1);
    J.block(dim, 3 * dim, dim, dim).diagonal().setConstant(-1);
    return J;
}

} // namespace ipc
