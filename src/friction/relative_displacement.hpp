#pragma once

#include <Eigen/Core>

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

template <
    typename DerivedDP0,
    typename DerivedDP1,
    typename T = typename DerivedDP0::Scalar>
inline MatrixMax<T, 3, 6> point_point_relative_displacement_jacobian(
    const Eigen::MatrixBase<DerivedDP0>& dp0,
    const Eigen::MatrixBase<DerivedDP1>& dp1)
{
    int n = dp0.size();
    assert(dp0.size() == n);
    MatrixMax<T, 3, 6> J(n, 2 * n);
    J.leftCols(n) = MatrixMax<T, 3, 6>::Identity(n, n);
    J.rightCols(n) = -MatrixMax<T, 3, 6>::Identity(n, n);
    return J;
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

template <
    typename DerivedDP,
    typename DerivedDE0,
    typename DerivedDE1,
    typename T>
inline MatrixMax<T, 3, 9> point_edge_relative_displacement_jacobian(
    const Eigen::MatrixBase<DerivedDP>& dp,
    const Eigen::MatrixBase<DerivedDE0>& de0,
    const Eigen::MatrixBase<DerivedDE1>& de1,
    const T& alpha)
{
    int n = dp.size();
    assert(de0.size() == n && de1.size() == n);
    MatrixMax<T, 3, 9> J = MatrixMax<T, 3, 9>::Zero(n, 3 * n);
    J.leftCols(n).diagonal().setOnes();
    J.middleCols(n, n).diagonal().setConstant(alpha - 1);
    J.rightCols(n).diagonal().setConstant(-alpha);
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

template <
    typename DerivedDEA0,
    typename DerivedDEA1,
    typename DerivedDEB0,
    typename DerivedDEB1,
    typename DerivedCoords,
    typename T = typename DerivedDEA0::Scalar>
inline MatrixMax<T, 3, 12> edge_edge_relative_displacement_jacobian(
    const Eigen::MatrixBase<DerivedDEA0>& dea0,
    const Eigen::MatrixBase<DerivedDEA1>& dea1,
    const Eigen::MatrixBase<DerivedDEB0>& deb0,
    const Eigen::MatrixBase<DerivedDEB1>& deb1,
    const Eigen::MatrixBase<DerivedCoords>& coords)
{
    int n = dea0.size();
    assert(dea1.size() == n && deb0.size() == n && deb1.size() == n);
    MatrixMax<T, 3, 12> J = MatrixMax<T, 3, 12>::Zero(n, 4 * n);
    J.leftCols(n).diagonal().setConstant(1 - coords[0]);
    J.middleCols(n, n).diagonal().setConstant(coords[0]);
    J.middleCols(2 * n, n).diagonal().setConstant(coords[1] - 1);
    J.rightCols(n).diagonal().setConstant(-coords[1]);
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

template <
    typename DerivedDP,
    typename DerivedDT0,
    typename DerivedDT1,
    typename DerivedDT2,
    typename DerivedCoords,
    typename T = typename DerivedDP::Scalar>
inline MatrixMax<T, 3, 12> point_triangle_relative_displacement_jacobian(
    const Eigen::MatrixBase<DerivedDP>& dp,
    const Eigen::MatrixBase<DerivedDT0>& dt0,
    const Eigen::MatrixBase<DerivedDT1>& dt1,
    const Eigen::MatrixBase<DerivedDT2>& dt2,
    const Eigen::MatrixBase<DerivedCoords>& coords)
{
    int n = dp.size();
    assert(dt0.size() == n && dt1.size() == n && dt2.size() == n);
    MatrixMax<T, 3, 12> J = MatrixMax<T, 3, 12>::Zero(n, 4 * n);
    J.leftCols(n).diagonal().setOnes();
    J.middleCols(n, n).diagonal().setConstant(coords[0] + coords[1] - 1);
    J.middleCols(2 * n, n).diagonal().setConstant(-coords[0]);
    J.rightCols(n).diagonal().setConstant(-coords[1]);
    return J;
}

} // namespace ipc
