#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <tuple>

namespace ipc {

// =============================================================================

/// @brief Computes the normalization and Jacobian of a vector.
/// @param x The input vector.
/// @return A tuple containing the normalized vector and its Jacobian.
std::tuple<VectorMax3d, MatrixMax3d>
normalization_and_jacobian(Eigen::ConstRef<VectorMax3d> x);

/// @brief Computes the Jacobian of the normalization operation.
/// @param x The input vector.
/// @return The Jacobian of the normalization operation.
inline MatrixMax3d normalization_jacobian(Eigen::ConstRef<VectorMax3d> x)
{
    return std::get<1>(normalization_and_jacobian(x));
}

/// @brief Computes the normalization, Jacobian, and Hessian of a vector.
/// @param t The input vector.
/// @return A tuple containing the normalized vector, its Jacobian, and its Hessian.
std::tuple<VectorMax3d, MatrixMax3d, std::array<MatrixMax3d, 3>>
normalization_and_jacobian_and_hessian(Eigen::ConstRef<VectorMax3d> x);

/// @brief Computes the Hessian of the normalization operation.
/// @param x The input vector.
/// @return The Hessian of the normalization operation.
inline std::array<MatrixMax3d, 3>
normalization_hessian(Eigen::ConstRef<VectorMax3d> x)
{
    return std::get<2>(normalization_and_jacobian_and_hessian(x));
}

// =============================================================================

/// @brief Cross product matrix for 3D vectors.
/// @param v Vector to create the cross product matrix for.
/// @return The cross product matrix of the vector.
inline Eigen::Matrix3d cross_product_matrix(Eigen::ConstRef<Eigen::Vector3d> v)
{
    Eigen::Matrix3d m;
    m << 0, -v(2), v(1), //
        v(2), 0, -v(0),  //
        -v(1), v(0), 0;
    return m;
}

/// @brief Computes the Jacobian of the cross product matrix.
/// @return The Jacobian of the cross product matrix.
Eigen::Matrix<double, 9, 3> cross_product_matrix_jacobian();

// =============================================================================

/**
 * \defgroup geometry Point-line normal
 * \brief Functions for computing a point-line normal and resp. Jacobians.
 * @{
 */

/// @brief Computes the unnormalized normal vector of a point-line pair.
/// @param p The vertex position.
/// @param e0 The start position of the line.
/// @param e1 The end position of the line.
/// @return The unnormalized normal vector.
VectorMax3d point_line_unnormalized_normal(
    Eigen::ConstRef<VectorMax3d> p,
    Eigen::ConstRef<VectorMax3d> e0,
    Eigen::ConstRef<VectorMax3d> e1);

/// @brief Computes the normal vector of a point-line pair.
/// @param p The vertex position.
/// @param e0 The start position of the line.
/// @param e1 The end position of the line.
/// @return The normal vector.
inline VectorMax3d point_line_normal(
    Eigen::ConstRef<VectorMax3d> p,
    Eigen::ConstRef<VectorMax3d> e0,
    Eigen::ConstRef<VectorMax3d> e1)
{
    return point_line_unnormalized_normal(p, e0, e1).normalized();
}

/// @brief Computes the Jacobian of the unnormalized normal vector of a point-line pair.
/// @param p The vertex position.
/// @param e0 The start position of the line.
/// @param e1 The end position of the line.
/// @return The Jacobian of the unnormalized normal vector.
MatrixMax<double, 3, 9> point_line_unnormalized_normal_jacobian(
    Eigen::ConstRef<VectorMax3d> p,
    Eigen::ConstRef<VectorMax3d> e0,
    Eigen::ConstRef<VectorMax3d> e1);

/// @brief Computes the Hessian of the unnormalized normal vector of a point-line pair.
/// @param p The vertex position.
/// @param e0 The start position of the line.
/// @param e1 The end position of the line.
/// @return The Hessian of the unnormalized normal vector of the point-line pair.
MatrixMax<double, 27, 9> point_line_unnormalized_normal_hessian(
    Eigen::ConstRef<VectorMax3d> p,
    Eigen::ConstRef<VectorMax3d> e0,
    Eigen::ConstRef<VectorMax3d> e1);

/// @brief Computes the Jacobian of the normal vector of a point-line pair.
/// @param p The vertex position.
/// @param e0 The start position of the line.
/// @param e1 The end position of the line.
/// @return The Jacobian of the normal vector.
inline MatrixMax<double, 3, 9> point_line_normal_jacobian(
    Eigen::ConstRef<VectorMax3d> p,
    Eigen::ConstRef<VectorMax3d> e0,
    Eigen::ConstRef<VectorMax3d> e1)
{
    // ∂n̂/∂x = ∂n̂/∂n * ∂n/∂x
    return normalization_jacobian(point_line_unnormalized_normal(p, e0, e1))
        * point_line_unnormalized_normal_jacobian(p, e0, e1);
}

/// @brief Computes the Hessian of the normal vector of a point-line pair.
/// @param p The vertex position.
/// @param e0 The start position of the line.
/// @param e1 The end position of the line.
/// @return The Hessian of the normal vector.
MatrixMax<double, 27, 9> point_line_normal_hessian(
    Eigen::ConstRef<VectorMax3d> p,
    Eigen::ConstRef<VectorMax3d> e0,
    Eigen::ConstRef<VectorMax3d> e1);

/** @} */

// =============================================================================

/**
 * \defgroup geometry Triangle normal
 * \brief Functions for computing a triangle's normal and resp. Jacobians.
 * @{
 */

/// @brief Computes the unnormalized normal vector of a triangle.
/// @param a The first vertex of the triangle.
/// @param b The second vertex of the triangle.
/// @param c The third vertex of the triangle.
/// @return The unnormalized normal vector of the triangle.
inline Eigen::Vector3d triangle_unnormalized_normal(
    Eigen::ConstRef<Eigen::Vector3d> a,
    Eigen::ConstRef<Eigen::Vector3d> b,
    Eigen::ConstRef<Eigen::Vector3d> c)
{
    return (b - a).cross(c - a);
}

/// @brief Computes the normal vector of a triangle.
/// @param a The first vertex of the triangle.
/// @param b The second vertex of the triangle.
/// @param c The third vertex of the triangle.
/// @return The normal vector of the triangle.
inline Eigen::Vector3d triangle_normal(
    Eigen::ConstRef<Eigen::Vector3d> a,
    Eigen::ConstRef<Eigen::Vector3d> b,
    Eigen::ConstRef<Eigen::Vector3d> c)
{
    return triangle_unnormalized_normal(a, b, c).normalized();
}

/// @brief Computes the Jacobian of the unnormalized normal vector of a triangle.
/// @param a The first vertex of the triangle.
/// @param b The second vertex of the triangle.
/// @param c The third vertex of the triangle.
/// @return The Jacobian of the unnormalized normal vector of the triangle.
inline Eigen::Matrix<double, 3, 9> triangle_unnormalized_normal_jacobian(
    Eigen::ConstRef<Eigen::Vector3d> a,
    Eigen::ConstRef<Eigen::Vector3d> b,
    Eigen::ConstRef<Eigen::Vector3d> c)
{
    Eigen::Matrix<double, 3, 9> J;
    J.middleCols<3>(0) = cross_product_matrix(c - b); // ∂n/∂a
    J.middleCols<3>(3) = cross_product_matrix(a - c); // ∂n/∂b
    J.middleCols<3>(6) = cross_product_matrix(b - a); // ∂n/∂c
    return J;
}

/// @brief Computes the Hessian of the unnormalized normal vector of a triangle.
/// @param a The first vertex of the triangle.
/// @param b The second vertex of the triangle.
/// @param c The third vertex of the triangle.
/// @return The Hessian of the unnormalized normal vector of the triangle.
Eigen::Matrix<double, 27, 9> triangle_unnormalized_normal_hessian(
    Eigen::ConstRef<Eigen::Vector3d> a,
    Eigen::ConstRef<Eigen::Vector3d> b,
    Eigen::ConstRef<Eigen::Vector3d> c);

/// @brief Computes the Jacobian of the normal vector of a triangle.
/// @param a The first vertex of the triangle.
/// @param b The second vertex of the triangle.
/// @param c The third vertex of the triangle.
/// @return The Jacobian of the normal vector of the triangle.
inline Eigen::Matrix<double, 3, 9> triangle_normal_jacobian(
    Eigen::ConstRef<Eigen::Vector3d> a,
    Eigen::ConstRef<Eigen::Vector3d> b,
    Eigen::ConstRef<Eigen::Vector3d> c)
{
    // ∂n̂/∂x = ∂n̂/∂n * ∂n/∂x
    return normalization_jacobian(triangle_unnormalized_normal(a, b, c))
        * triangle_unnormalized_normal_jacobian(a, b, c);
}

/// @brief Computes the Hessian of the normal vector of a triangle.
/// @param a The first vertex of the triangle.
/// @param b The second vertex of the triangle.
/// @param c The third vertex of the triangle.
/// @return The Hessian of the normal vector of the triangle.
Eigen::Matrix<double, 27, 9> triangle_normal_hessian(
    Eigen::ConstRef<Eigen::Vector3d> a,
    Eigen::ConstRef<Eigen::Vector3d> b,
    Eigen::ConstRef<Eigen::Vector3d> c);

/** @} */

// =============================================================================

/**
 * \defgroup geometry Line-line normal
 * \brief Functions for computing a line-line normal and resp. Jacobians.
 * @{
 */

/// @brief Computes the unnormalized normal vector of two lines.
/// @param ea0 The first vertex of the first line.
/// @param ea1 The second vertex of the first line.
/// @param eb0 The first vertex of the second line.
/// @param eb1 The second vertex of the second line.
/// @return The unnormalized normal vector of the two lines.
inline Eigen::Vector3d line_line_unnormalized_normal(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
    return (ea1 - ea0).cross(eb1 - eb0);
}

/// @brief Computes the normal vector of two lines.
/// @param ea0 The first vertex of the first line.
/// @param ea1 The second vertex of the first line.
/// @param eb0 The first vertex of the second line.
/// @param eb1 The second vertex of the second line.
/// @return The normal vector of the two lines.
inline Eigen::Vector3d line_line_normal(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
    return line_line_unnormalized_normal(ea0, ea1, eb0, eb1).normalized();
}

/// @brief Computes the Jacobian of the unnormalized normal vector of two lines.
/// @param ea0 The first vertex of the first line.
/// @param ea1 The second vertex of the first line.
/// @param eb0 The first vertex of the second line.
/// @param eb1 The second vertex of the second line.
/// @return The Jacobian of the unnormalized normal vector of the two lines.
inline Eigen::Matrix<double, 3, 12> line_line_unnormalized_normal_jacobian(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
    Eigen::Matrix<double, 3, 12> J;
    J.middleCols<3>(0) = cross_product_matrix(eb1 - eb0);
    J.middleCols<3>(3) = cross_product_matrix(eb0 - eb1);
    J.middleCols<3>(6) = cross_product_matrix(ea0 - ea1);
    J.middleCols<3>(9) = cross_product_matrix(ea1 - ea0);
    return J;
}

/// @brief Computes the Jacobian of the normal vector of two lines.
/// @param ea0 The first vertex of the first line.
/// @param ea1 The second vertex of the first line.
/// @param eb0 The first vertex of the second line.
/// @param eb1 The second vertex of the second line.
/// @return The Jacobian of the normal vector of the two lines.
inline Eigen::Matrix<double, 3, 12> line_line_normal_jacobian(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
    // ∂n̂/∂x = ∂n̂/∂n * ∂n/∂x
    return normalization_jacobian(
               line_line_unnormalized_normal(ea0, ea1, eb0, eb1))
        * line_line_unnormalized_normal_jacobian(ea0, ea1, eb0, eb1);
}

/// @brief Computes the Hessian of the unnormalized normal vector of two lines.
/// @param ea0 The first vertex of the first line.
/// @param ea1 The second vertex of the first line.
/// @param eb0 The first vertex of the second line.
/// @param eb1 The second vertex of the second line.
/// @return The Hessian of the unnormalized normal vector of the two lines.
Eigen::Matrix<double, 36, 12> line_line_unnormalized_normal_hessian(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1);

/// @brief Computes the Hessian of the normal vector of two lines.
/// @param ea0 The first vertex of the first line.
/// @param ea1 The second vertex of the first line.
/// @param eb0 The first vertex of the second line.
/// @param eb1 The second vertex of the second line.
/// @return The Hessian of the normal vector of the two lines.
Eigen::Matrix<double, 36, 12> line_line_normal_hessian(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1);

/** @} */

} // namespace ipc