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

VectorMax<double, 18> point_point_relative_velocity_dx_dbeta(const int dim)
{
    // Γ is constant (does not depend on β), so the derivative is zero.
    return VectorMax<double, 18>::Zero(2 * dim * dim);
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

// Γ(α) = [I, (α-1)I, -αI]  (dim × 3·dim)
// ∂Γ/∂α = [0, I, -I]
//
// Stored as vec(∂Γ/∂α) in column-major order (3rd-order convention).
// Result is a column vector of size dim × 3·dim = 3·dim².
// For a (dim × ndof) matrix M, element M(r,c) maps to vec index c·dim + r.
VectorMax<double, 27>
point_edge_relative_velocity_dx_dbeta(const int dim, const double alpha)
{
    const int ndof = 3 * dim;
    VectorMax<double, 27> J = VectorMax<double, 27>::Zero(dim * ndof);
    for (int i = 0; i < dim; ++i) {
        // I block at cols [dim, 2·dim)
        J[(dim + i) * dim + i] = 1;
        // -I block at cols [2·dim, 3·dim)
        J[(2 * dim + i) * dim + i] = -1;
    }
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

// Γ(β₁,β₂) = [(1-β₁)I, β₁I, (β₂-1)I, -β₂I]  (3 × 12)
// ∂Γ/∂β₁ = [-I, I, 0, 0]
// ∂Γ/∂β₂ = [ 0, 0, I,-I]
//
// Stored as [vec(∂Γ/∂β₁) | vec(∂Γ/∂β₂)] in column-major order (3rd-order
// convention). Result shape: (36, 2). For a (3 × 12) matrix M, element M(r,c)
// maps to vec index c·3 + r.
Eigen::Matrix<double, 36, 2>
edge_edge_relative_velocity_dx_dbeta(Eigen::ConstRef<Eigen::Vector2d> coords)
{
    constexpr int dim = 3;
    Eigen::Matrix<double, 36, 2> J = Eigen::Matrix<double, 36, 2>::Zero();
    for (int i = 0; i < dim; ++i) {
        // wrt β₁: -I at cols [0,3), I at cols [3,6)
        J((0 + i) * dim + i, 0) = -1;
        J((3 + i) * dim + i, 0) = 1;
        // wrt β₂: I at cols [6,9), -I at cols [9,12)
        J((6 + i) * dim + i, 1) = 1;
        J((9 + i) * dim + i, 1) = -1;
    }
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

// Γ(β₁,β₂) = [I, (β₁+β₂-1)I, -β₁I, -β₂I]  (3 × 12)
// ∂Γ/∂β₁ = [0, I, -I, 0]
// ∂Γ/∂β₂ = [0, I,  0,-I]
//
// Stored as [vec(∂Γ/∂β₁) | vec(∂Γ/∂β₂)] in column-major order (3rd-order
// convention). Result shape: (36, 2). For a (3 × 12) matrix M, element M(r,c)
// maps to vec index c·3 + r.
Eigen::Matrix<double, 36, 2> point_triangle_relative_velocity_dx_dbeta(
    Eigen::ConstRef<Eigen::Vector2d> coords)
{
    constexpr int dim = 3;
    Eigen::Matrix<double, 36, 2> J = Eigen::Matrix<double, 36, 2>::Zero();
    for (int i = 0; i < dim; ++i) {
        // wrt β₁: I at cols [3,6), -I at cols [6,9)
        J((3 + i) * dim + i, 0) = 1;
        J((6 + i) * dim + i, 0) = -1;
        // wrt β₂: I at cols [3,6), -I at cols [9,12)
        J((3 + i) * dim + i, 1) = 1;
        J((9 + i) * dim + i, 1) = -1;
    }
    return J;
}

} // namespace ipc