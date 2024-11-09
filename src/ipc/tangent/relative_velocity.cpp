#include "relative_velocity.hpp"

namespace ipc {

// ============================================================================
// Point - Point

VectorMax3d point_point_relative_velocity(
    const Eigen::Ref<const VectorMax3d>& dp0,
    const Eigen::Ref<const VectorMax3d>& dp1)
{
    return dp0 - dp1;
}

MatrixMax<double, 3, 6> point_point_relative_velocity_matrix(const int dim)
{
    MatrixMax<double, 3, 6> J(dim, 2 * dim);
    J.leftCols(dim) = MatrixMax3d::Identity(dim, dim);
    J.rightCols(dim) = -MatrixMax3d::Identity(dim, dim);
    return J;
}

MatrixMax<double, 3, 6>
point_point_relative_velocity_matrix_jacobian(const int dim)
{
    return MatrixMax<double, 3, 6>::Zero(dim, 2 * dim);
}

// ============================================================================
// Point - Edge

VectorMax3d point_edge_relative_velocity(
    const Eigen::Ref<const VectorMax3d>& dp,
    const Eigen::Ref<const VectorMax3d>& de0,
    const Eigen::Ref<const VectorMax3d>& de1,
    const double alpha)
{
    return dp - ((de1 - de0) * alpha + de0);
}

MatrixMax<double, 3, 9>
point_edge_relative_velocity_matrix(const int dim, const double alpha)
{
    MatrixMax<double, 3, 9> J = MatrixMax<double, 3, 9>::Zero(dim, 3 * dim);
    J.leftCols(dim).diagonal().setOnes();
    J.middleCols(dim, dim).diagonal().setConstant(alpha - 1);
    J.rightCols(dim).diagonal().setConstant(-alpha);
    return J;
}

MatrixMax<double, 3, 9>
point_edge_relative_velocity_matrix_jacobian(const int dim, const double alpha)
{
    MatrixMax<double, 3, 9> J = MatrixMax<double, 3, 9>::Zero(dim, 3 * dim);
    J.middleCols(dim, dim).diagonal().setConstant(1);
    J.rightCols(dim).diagonal().setConstant(-1);
    return J;
}

// ============================================================================
// Edge - Edge

Eigen::Vector3d edge_edge_relative_velocity(
    const Eigen::Ref<const Eigen::Vector3d>& dea0,
    const Eigen::Ref<const Eigen::Vector3d>& dea1,
    const Eigen::Ref<const Eigen::Vector3d>& deb0,
    const Eigen::Ref<const Eigen::Vector3d>& deb1,
    const Eigen::Ref<const Eigen::Vector2d>& coords)
{
    // closest_point_a_velocity - closest_point_b_velocity
    return ((dea1 - dea0) * coords[0] + dea0)
        - ((deb1 - deb0) * coords[1] + deb0);
}

MatrixMax<double, 3, 12> edge_edge_relative_velocity_matrix(
    const int dim, const Eigen::Ref<const Eigen::Vector2d>& coords)
{
    MatrixMax<double, 3, 12> J = MatrixMax<double, 3, 12>::Zero(dim, 4 * dim);
    J.leftCols(dim).diagonal().setConstant(1 - coords[0]);
    J.middleCols(dim, dim).diagonal().setConstant(coords[0]);
    J.middleCols(2 * dim, dim).diagonal().setConstant(coords[1] - 1);
    J.rightCols(dim).diagonal().setConstant(-coords[1]);
    return J;
}

MatrixMax<double, 6, 12> edge_edge_relative_velocity_matrix_jacobian(
    const int dim, const Eigen::Ref<const Eigen::Vector2d>& coords)
{
    MatrixMax<double, 6, 12> J =
        MatrixMax<double, 6, 12>::Zero(2 * dim, 4 * dim);
    // wrt β₁
    J.block(0, 0, dim, dim).diagonal().setConstant(-1);
    J.block(0, dim, dim, dim).diagonal().setConstant(1);
    // wrt β₂
    J.block(dim, 2 * dim, dim, dim).diagonal().setConstant(1);
    J.block(dim, 3 * dim, dim, dim).diagonal().setConstant(-1);
    return J;
}

// ============================================================================
// Point - Triangle

Eigen::Vector3d point_triangle_relative_velocity(
    const Eigen::Ref<const Eigen::Vector3d>& dp,
    const Eigen::Ref<const Eigen::Vector3d>& dt0,
    const Eigen::Ref<const Eigen::Vector3d>& dt1,
    const Eigen::Ref<const Eigen::Vector3d>& dt2,
    const Eigen::Ref<const Eigen::Vector2d>& coords)
{
    // Compute the velocity of the closest point and subtract it from the
    // points velocity.
    return dp - (dt0 + coords[0] * (dt1 - dt0) + coords[1] * (dt2 - dt0));
}

MatrixMax<double, 3, 12> point_triangle_relative_velocity_matrix(
    const int dim, const Eigen::Ref<const Eigen::Vector2d>& coords)
{
    MatrixMax<double, 3, 12> J = MatrixMax<double, 3, 12>::Zero(dim, 4 * dim);
    J.leftCols(dim).diagonal().setOnes();
    J.middleCols(dim, dim).diagonal().setConstant(coords[0] + coords[1] - 1);
    J.middleCols(2 * dim, dim).diagonal().setConstant(-coords[0]);
    J.rightCols(dim).diagonal().setConstant(-coords[1]);
    return J;
}

MatrixMax<double, 6, 12> point_triangle_relative_velocity_matrix_jacobian(
    const int dim, const Eigen::Ref<const Eigen::Vector2d>& coords)
{
    MatrixMax<double, 6, 12> J =
        MatrixMax<double, 6, 12>::Zero(2 * dim, 4 * dim);
    // wrt β₁
    J.block(0, dim, dim, dim).diagonal().setConstant(1);
    J.block(0, 2 * dim, dim, dim).diagonal().setConstant(-1);
    // wrt β₂
    J.block(dim, dim, dim, dim).diagonal().setConstant(1);
    J.block(dim, 3 * dim, dim, dim).diagonal().setConstant(-1);
    return J;
}

} // namespace ipc
