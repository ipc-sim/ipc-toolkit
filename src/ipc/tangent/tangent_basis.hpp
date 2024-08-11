#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

// ============================================================================
// Point - Point

/// @brief Compute a basis for the space tangent to the point-point pair.
/// @param p0 First point
/// @param p1 Second point
/// @return A 3x2 matrix whose columns are the basis vectors.
MatrixMax<double, 3, 2> point_point_tangent_basis(
    const Eigen::Ref<const VectorMax3d>& p0,
    const Eigen::Ref<const VectorMax3d>& p1);

/// @brief Compute the Jacobian of the tangent basis for the point-point pair.
/// @param p0 First point
/// @param p1 Second point
/// @return A 6*3x2 matrix whose columns are the basis vectors.
MatrixMax<double, 18, 2> point_point_tangent_basis_jacobian(
    const Eigen::Ref<const VectorMax3d>& p0,
    const Eigen::Ref<const VectorMax3d>& p1);

// ============================================================================
// Point - Edge

/// @brief Compute a basis for the space tangent to the point-edge pair.
/// @param p Point
/// @param e0 First edge point
/// @param e1 Second edge point
/// @return A 3x2 matrix whose columns are the basis vectors.
MatrixMax<double, 3, 2> point_edge_tangent_basis(
    const Eigen::Ref<const VectorMax3d>& p,
    const Eigen::Ref<const VectorMax3d>& e0,
    const Eigen::Ref<const VectorMax3d>& e1);

/// @brief Compute the Jacobian of the tangent basis for the point-edge pair.
/// @param p Point
/// @param e0 First edge point
/// @param e1 Second edge point
/// @return A 9*3x2 matrix whose columns are the basis vectors.
MatrixMax<double, 27, 2> point_edge_tangent_basis_jacobian(
    const Eigen::Ref<const VectorMax3d>& p,
    const Eigen::Ref<const VectorMax3d>& e0,
    const Eigen::Ref<const VectorMax3d>& e1);

// ============================================================================
// Edge - Edge

/// @brief Compute a basis for the space tangent to the edge-edge pair.
/// @param ea0 First point of the first edge
/// @param ea1 Second point of the first edge
/// @param eb0 First point of the second edge
/// @param eb1 Second point of the second edge
/// @return A 3x2 matrix whose columns are the basis vectors.
Eigen::Matrix<double, 3, 2> edge_edge_tangent_basis(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1);

/// @brief Compute the Jacobian of the tangent basis for the edge-edge pair.
/// @param ea0 First point of the first edge
/// @param ea1 Second point of the first edge
/// @param eb0 First point of the second edge
/// @param eb1 Second point of the second edge
/// @return A 12*3x2 matrix whose columns are the basis vectors.
Eigen::Matrix<double, 36, 2> edge_edge_tangent_basis_jacobian(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1);

// ============================================================================
// Point - Triangle

/// @brief Compute a basis for the space tangent to the point-triangle pair.
///
/// \f\[
///     \begin{bmatrix}
///     \frac{t_1 - t_0}{\|t_1 - t_0\|} & \frac{((t_1 - t_0)\times(t_2 - t_0))
///     \times(t_1 - t_0)}{\|((t_1 - t_0)\times(t_2 - t_0))\times(t_1 - t_0)\|}
///     \end{bmatrix}
/// \f\]
///
/// @param p Point
/// @param t0 Triangle's first vertex
/// @param t1 Triangle's second vertex
/// @param t2 Triangle's third vertex
/// @return A 3x2 matrix whose columns are the basis vectors.
Eigen::Matrix<double, 3, 2> point_triangle_tangent_basis(
    const Eigen::Ref<const Eigen::Vector3d>& p,
    const Eigen::Ref<const Eigen::Vector3d>& t0,
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2);

/// @brief Compute the Jacobian of the tangent basis for the point-triangle pair.
/// @param p Point
/// @param t0 Triangle's first vertex
/// @param t1 Triangle's second vertex
/// @param t2 Triangle's third vertex
/// @return A 12*3x2 matrix whose columns are the basis vectors.
Eigen::Matrix<double, 36, 2> point_triangle_tangent_basis_jacobian(
    const Eigen::Ref<const Eigen::Vector3d>& p,
    const Eigen::Ref<const Eigen::Vector3d>& t0,
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2);

// ============================================================================

namespace autogen {
    // J is (8×1) flattened in column-major order
    void point_point_tangent_basis_2D_jacobian(
        double p0_x, double p0_y, double p1_x, double p1_y, double J[8]);

    // J is (18×2) flattened in column-major order
    void point_point_tangent_basis_3D_jacobian(
        double p0_x,
        double p0_y,
        double p0_z,
        double p1_x,
        double p1_y,
        double p1_z,
        double J[36]);

    // J is (12×1) flattened in column-major order
    void point_edge_tangent_basis_2D_jacobian(
        double p_x,
        double p_y,
        double e0_x,
        double e0_y,
        double e1_x,
        double e1_y,
        double J[12]);

    // J is (27×2) flattened in column-major order
    void point_edge_tangent_basis_3D_jacobian(
        double p_x,
        double p_y,
        double p_z,
        double e0_x,
        double e0_y,
        double e0_z,
        double e1_x,
        double e1_y,
        double e1_z,
        double J[54]);

    // J is (36×2) flattened in column-major order
    void edge_edge_tangent_basis_jacobian(
        double ea0_x,
        double ea0_y,
        double ea0_z,
        double ea1_x,
        double ea1_y,
        double ea1_z,
        double eb0_x,
        double eb0_y,
        double eb0_z,
        double eb1_x,
        double eb1_y,
        double eb1_z,
        double J[72]);

    // J is (36×2) flattened in column-major order
    void point_triangle_tangent_basis_jacobian(
        double p_x,
        double p_y,
        double p_z,
        double t0_x,
        double t0_y,
        double t0_z,
        double t1_x,
        double t1_y,
        double t1_z,
        double t2_x,
        double t2_y,
        double t2_z,
        double J[72]);
} // namespace autogen

} // namespace ipc
