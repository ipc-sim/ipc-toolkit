#include "relative_velocity.hpp"

namespace ipc {

// ============================================================================
// Point - Point

VectorMax3d point_point_relative_velocity(
    Eigen::ConstRef<VectorMax3d> dp0, Eigen::ConstRef<VectorMax3d> dp1)
{
    return dp0 - dp1;
}

MatrixMax<double, 3, 6> point_point_relative_velocity_jacobian(const int dim)
{
    MatrixMax<double, 3, 6> J(dim, 2 * dim);
    J.leftCols(dim) = MatrixMax3d::Identity(dim, dim);
    J.rightCols(dim) = -MatrixMax3d::Identity(dim, dim);
    return J;
}

MatrixMax<double, 3, 6> point_point_relative_velocity_dx_dbeta(const int dim)
{
    return MatrixMax<double, 3, 6>::Zero(dim, 2 * dim);
}

// ============================================================================
// Point - Edge

VectorMax3d point_edge_relative_velocity(
    Eigen::ConstRef<VectorMax3d> dp,
    Eigen::ConstRef<VectorMax3d> de0,
    Eigen::ConstRef<VectorMax3d> de1,
    const double alpha)
{
    return dp - ((de1 - de0) * alpha + de0);
}

MatrixMax<double, 3, 9>
point_edge_relative_velocity_jacobian(const int dim, const double alpha)
{
    MatrixMax<double, 3, 9> J = MatrixMax<double, 3, 9>::Zero(dim, 3 * dim);
    J.leftCols(dim).diagonal().setOnes();
    J.middleCols(dim, dim).diagonal().setConstant(alpha - 1);
    J.rightCols(dim).diagonal().setConstant(-alpha);
    return J;
}

MatrixMax<double, 3, 9>
point_edge_relative_velocity_dx_dbeta(const int dim, const double alpha)
{
    MatrixMax<double, 3, 9> J = MatrixMax<double, 3, 9>::Zero(dim, 3 * dim);
    J.middleCols(dim, dim).diagonal().setConstant(1);
    J.rightCols(dim).diagonal().setConstant(-1);
    return J;
}

// ============================================================================
// Edge - Edge

Eigen::Vector3d edge_edge_relative_velocity(
    Eigen::ConstRef<Eigen::Vector3d> dea0,
    Eigen::ConstRef<Eigen::Vector3d> dea1,
    Eigen::ConstRef<Eigen::Vector3d> deb0,
    Eigen::ConstRef<Eigen::Vector3d> deb1,
    Eigen::ConstRef<Eigen::Vector2d> coords)
{
    // closest_point_a_velocity - closest_point_b_velocity
    return ((dea1 - dea0) * coords[0] + dea0)
        - ((deb1 - deb0) * coords[1] + deb0);
}

Eigen::Matrix<double, 3, 12>
edge_edge_relative_velocity_jacobian(Eigen::ConstRef<Eigen::Vector2d> coords)
{
    Eigen::Matrix<double, 3, 12> J = Eigen::Matrix<double, 3, 12>::Zero();
    J.leftCols<3>().diagonal().setConstant(1 - coords[0]);
    J.middleCols<3>(3).diagonal().setConstant(coords[0]);
    J.middleCols<3>(6).diagonal().setConstant(coords[1] - 1);
    J.rightCols<3>().diagonal().setConstant(-coords[1]);
    return J;
}

Eigen::Matrix<double, 3, 24>
edge_edge_relative_velocity_dx_dbeta(Eigen::ConstRef<Eigen::Vector2d> coords)
{
    Eigen::Matrix<double, 3, 24> J = Eigen::Matrix<double, 3, 24>::Zero();
    // wrt β₁
    J.middleCols<3>(0).diagonal().setConstant(-1);
    J.middleCols<3>(3).diagonal().setConstant(1);
    // // wrt β₂
    J.middleCols<3>(18).diagonal().setConstant(1);
    J.middleCols<3>(21).diagonal().setConstant(-1);
    return J;
}

// ============================================================================
// Point - Triangle

Eigen::Vector3d point_triangle_relative_velocity(
    Eigen::ConstRef<Eigen::Vector3d> dp,
    Eigen::ConstRef<Eigen::Vector3d> dt0,
    Eigen::ConstRef<Eigen::Vector3d> dt1,
    Eigen::ConstRef<Eigen::Vector3d> dt2,
    Eigen::ConstRef<Eigen::Vector2d> coords)
{
    // Compute the velocity of the closest point and subtract it from the
    // points velocity.
    return dp - (dt0 + coords[0] * (dt1 - dt0) + coords[1] * (dt2 - dt0));
}

Eigen::Matrix<double, 3, 12> point_triangle_relative_velocity_jacobian(
    Eigen::ConstRef<Eigen::Vector2d> coords)
{
    Eigen::Matrix<double, 3, 12> J = Eigen::Matrix<double, 3, 12>::Zero();
    J.leftCols<3>().diagonal().setOnes();
    J.middleCols<3>(3).diagonal().setConstant(coords[0] + coords[1] - 1);
    J.middleCols<3>(6).diagonal().setConstant(-coords[0]);
    J.rightCols<3>().diagonal().setConstant(-coords[1]);
    return J;
}

Eigen::Matrix<double, 3, 24> point_triangle_relative_velocity_dx_dbeta(
    Eigen::ConstRef<Eigen::Vector2d> coords)
{
    Eigen::Matrix<double, 3, 24> J = Eigen::Matrix<double, 3, 24>::Zero();
    // wrt β₁
    J.middleCols<3>(3).diagonal().setConstant(1);
    J.middleCols<3>(6).diagonal().setConstant(-1);
    // wrt β₂
    J.middleCols<3>(15).diagonal().setConstant(1);
    J.middleCols<3>(21).diagonal().setConstant(-1);

    return J;
}

} // namespace ipc
