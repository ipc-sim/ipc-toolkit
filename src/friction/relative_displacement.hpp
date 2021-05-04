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
    typename DerivedDisp,
    typename DerivedBasis,
    typename T = typename DerivedDisp::Scalar>
inline VectorMax6<T> point_point_relative_mesh_displacements(
    const Eigen::MatrixBase<DerivedDisp>& tangent_relative_displacement,
    const Eigen::MatrixBase<DerivedBasis>& basis)
{
    int dim = basis.rows();
    VectorMax6<T> mesh_displacements(2 * dim);
    mesh_displacements.head(dim) = basis * tangent_relative_displacement;
    mesh_displacements.tail(dim) = -mesh_displacements.head(dim);
    return mesh_displacements;
}

template <typename DerivedBasis>
inline MatrixMax<typename DerivedBasis::Scalar, 2, 6>
point_point_TT(const Eigen::MatrixBase<DerivedBasis>& basis)
{
    MatrixMax<typename DerivedBasis::Scalar, 2, 6> TT(
        basis.cols(), 2 * basis.rows());
    TT.leftCols(basis.rows()) = basis.transpose();
    TT.rightCols(basis.rows()) = -basis.transpose();
    return TT;
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

template <typename DerivedDisp, typename DerivedBasis, typename T>
inline VectorMax9<T> point_edge_relative_mesh_displacements(
    const Eigen::MatrixBase<DerivedDisp>& tangent_relative_displacement,
    const Eigen::MatrixBase<DerivedBasis>& basis,
    const T& alpha)
{
    int dim = basis.rows();
    VectorMax9<T> mesh_displacements(3 * dim);
    mesh_displacements.head(dim) = basis * tangent_relative_displacement;
    mesh_displacements.segment(1 * dim, dim) =
        (alpha - 1.0) * mesh_displacements.head(dim);
    mesh_displacements.tail(dim) = -alpha * mesh_displacements.head(dim);
    return mesh_displacements;
}

template <typename DerivedBasis, typename T>
inline MatrixMax<typename DerivedBasis::Scalar, 2, 9>
point_edge_TT(const Eigen::MatrixBase<DerivedBasis>& basis, const T& alpha)
{
    MatrixMax<typename DerivedBasis::Scalar, 2, 9> TT(
        basis.cols(), 3 * basis.rows());
    TT.leftCols(basis.rows()) = basis.transpose();
    TT.middleCols(basis.rows(), basis.rows()) =
        (alpha - 1.0) * basis.transpose();
    TT.rightCols(basis.rows()) = -alpha * basis.transpose();
    return TT;
}

///////////////////////////////////////////////////////////////////////////////
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
inline VectorMax12<T> edge_edge_relative_mesh_displacements(
    const Eigen::MatrixBase<DerivedDisp>& tangent_relative_displacement,
    const Eigen::MatrixBase<DerivedBasis>& basis,
    const Eigen::MatrixBase<DerivedBeta>& gamma)
{
    VectorMax3<T> rel_disp = basis * tangent_relative_displacement;
    int dim = rel_disp.size();
    VectorMax12<T> mesh_displacements(4 * dim);
    mesh_displacements.head(dim) = (1.0 - gamma[0]) * rel_disp; // dea0
    mesh_displacements.segment(dim, dim) = gamma[0] * rel_disp; // dea1
    mesh_displacements.segment(2 * dim, dim) =
        (gamma[1] - 1.0) * rel_disp;                     // deb0
    mesh_displacements.tail(dim) = -gamma[1] * rel_disp; // deb1
    return mesh_displacements;
}

template <typename DerivedBasis, typename DerivedGamma>
inline MatrixMax<typename DerivedBasis::Scalar, 2, 12> edge_edge_TT(
    const Eigen::MatrixBase<DerivedBasis>& basis,
    const Eigen::MatrixBase<DerivedGamma>& gamma)
{
    MatrixMax<typename DerivedBasis::Scalar, 2, 12> TT(
        basis.cols(), 4 * basis.rows());
    TT.leftCols(basis.rows()) = (1.0 - gamma[0]) * basis.transpose();
    TT.middleCols(basis.rows(), basis.rows()) = gamma[0] * basis.transpose();
    TT.middleCols(2 * basis.rows(), basis.rows()) =
        (gamma[1] - 1.0) * basis.transpose();
    TT.rightCols(basis.rows()) = -gamma[1] * basis.transpose();
    return TT;
}

///////////////////////////////////////////////////////////////////////////////
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
inline VectorMax12<T> point_triangle_relative_mesh_displacements(
    const Eigen::MatrixBase<DerivedDisp>& tangent_relative_displacement,
    const Eigen::MatrixBase<DerivedBasis>& basis,
    const Eigen::MatrixBase<DerivedBeta>& beta)
{
    int dim = basis.rows();
    VectorMax12<T> mesh_displacements(4 * dim);
    mesh_displacements.head(dim) = basis * tangent_relative_displacement;
    mesh_displacements.segment(dim, dim) =
        (-1 + beta[0] + beta[1]) * mesh_displacements.head(dim);
    mesh_displacements.segment(2 * dim, dim) =
        -beta[0] * mesh_displacements.head(dim);
    mesh_displacements.tail(dim) = -beta[1] * mesh_displacements.head(dim);
    return mesh_displacements;
}

template <typename DerivedBasis, typename DerivedBeta>
inline MatrixMax<typename DerivedBasis::Scalar, 2, 12> point_triangle_TT(
    const Eigen::MatrixBase<DerivedBasis>& basis,
    const Eigen::MatrixBase<DerivedBeta>& beta)
{
    MatrixMax<typename DerivedBasis::Scalar, 2, 12> TT(
        basis.cols(), 4 * basis.rows());
    TT.leftCols(basis.rows()) = basis.transpose();
    TT.middleCols(basis.rows(), basis.rows()) =
        (-1 + beta[0] + beta[1]) * basis.transpose();
    TT.middleCols(2 * basis.rows(), basis.rows()) =
        -beta[0] * basis.transpose();
    TT.rightCols(basis.rows()) = -beta[1] * basis.transpose();
    return TT;
}

} // namespace ipc
