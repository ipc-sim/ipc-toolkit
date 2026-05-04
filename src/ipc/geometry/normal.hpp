#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <tuple>

namespace ipc {

// =============================================================================

/// @brief Computes the normalization and Jacobian of a vector.
/// @param x The input vector.
/// @return A tuple containing the normalized vector and its Jacobian.
template <typename T>
std::tuple<VectorMax3<T>, MatrixMax3<T>>
normalization_and_jacobian(Eigen::ConstRef<VectorMax3<T>> x);

/// @brief Computes the Jacobian of the normalization operation.
/// @param x The input vector.
/// @return The Jacobian of the normalization operation.
template <typename T>
inline MatrixMax3<T> normalization_jacobian(Eigen::ConstRef<VectorMax3<T>> x)
{
    return std::get<1>(normalization_and_jacobian(x));
}

/// @brief Computes the normalization, Jacobian, and Hessian of a vector.
/// @param x The input vector.
/// @return A tuple of the normalized vector, its Jacobian, and its Hessian.
template <typename T>
std::tuple<VectorMax3<T>, MatrixMax3<T>, std::array<MatrixMax3<T>, 3>>
normalization_and_jacobian_and_hessian(Eigen::ConstRef<VectorMax3<T>> x);

/// @brief Computes the Hessian of the normalization operation.
/// @param x The input vector.
/// @return The Hessian of the normalization operation.
template <typename T>
inline std::array<MatrixMax3<T>, 3>
normalization_hessian(Eigen::ConstRef<VectorMax3<T>> x)
{
    return std::get<2>(normalization_and_jacobian_and_hessian(x));
}

// =============================================================================

/// @brief Cross product matrix for 3D vectors.
/// @param v Vector to create the cross product matrix for.
/// @return The cross product matrix of the vector.
template <typename T>
inline Eigen::Matrix3<T>
cross_product_matrix(Eigen::ConstRef<Eigen::Vector3<T>> v)
{
    Eigen::Matrix3<T> m;
    // clang-format off
    m << T(0), -v(2), v(1),
         v(2),  T(0), -v(0),
        -v(1),  v(0),  T(0);
    // clang-format on
    return m;
}

/// @brief Computes the Jacobian of the cross product matrix.
/// @return The Jacobian of the cross product matrix.
template <typename T> Eigen::Matrix<T, 9, 3> cross_product_matrix_jacobian();

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
template <typename T>
VectorMax3<T> point_line_unnormalized_normal(
    Eigen::ConstRef<VectorMax3<T>> p,
    Eigen::ConstRef<VectorMax3<T>> e0,
    Eigen::ConstRef<VectorMax3<T>> e1);

/// @brief Computes the normal vector of a point-line pair.
/// @param p The vertex position.
/// @param e0 The start position of the line.
/// @param e1 The end position of the line.
/// @return The normal vector.
template <typename T>
inline VectorMax3<T> point_line_normal(
    Eigen::ConstRef<VectorMax3<T>> p,
    Eigen::ConstRef<VectorMax3<T>> e0,
    Eigen::ConstRef<VectorMax3<T>> e1)
{
    return point_line_unnormalized_normal(p, e0, e1).normalized();
}

/// @brief Computes the Jacobian of the unnormalized normal vector of a
/// point-line pair.
/// @param p The vertex position.
/// @param e0 The start position of the line.
/// @param e1 The end position of the line.
/// @return The Jacobian of the unnormalized normal vector.
template <typename T>
MatrixMax<T, 3, 9> point_line_unnormalized_normal_jacobian(
    Eigen::ConstRef<VectorMax3<T>> p,
    Eigen::ConstRef<VectorMax3<T>> e0,
    Eigen::ConstRef<VectorMax3<T>> e1);
/// @brief Computes the Hessian of the unnormalized normal vector of a
/// point-line pair.
/// @param p The vertex position.
/// @param e0 The start position of the line.
/// @param e1 The end position of the line.
/// @return The Hessian of the unnormalized normal vector of the point-line
/// pair.
template <typename T>
MatrixMax<T, 27, 9> point_line_unnormalized_normal_hessian(
    Eigen::ConstRef<VectorMax3<T>> p,
    Eigen::ConstRef<VectorMax3<T>> e0,
    Eigen::ConstRef<VectorMax3<T>> e1);

/// @brief Computes the Jacobian of the normal vector of a point-line pair.
/// @param p The vertex position.
/// @param e0 The start position of the line.
/// @param e1 The end position of the line.
/// @return The Jacobian of the normal vector.
template <typename T>
inline MatrixMax<T, 3, 9> point_line_normal_jacobian(
    Eigen::ConstRef<VectorMax3<T>> p,
    Eigen::ConstRef<VectorMax3<T>> e0,
    Eigen::ConstRef<VectorMax3<T>> e1)
{
    // ∂n̂/∂x = ∂n̂/∂n * ∂n/∂x
    return normalization_jacobian<T>(point_line_unnormalized_normal(p, e0, e1))
        * point_line_unnormalized_normal_jacobian(p, e0, e1);
}

/// @brief Computes the Hessian of the normal vector of a point-line pair.
/// @param p The vertex position.
/// @param e0 The start position of the line.
/// @param e1 The end position of the line.
/// @return The Hessian of the normal vector.
template <typename T>
MatrixMax<T, 27, 9> point_line_normal_hessian(
    Eigen::ConstRef<VectorMax3<T>> p,
    Eigen::ConstRef<VectorMax3<T>> e0,
    Eigen::ConstRef<VectorMax3<T>> e1);

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
template <typename T>
inline Eigen::Vector3<T> triangle_unnormalized_normal(
    Eigen::ConstRef<Eigen::Vector3<T>> a,
    Eigen::ConstRef<Eigen::Vector3<T>> b,
    Eigen::ConstRef<Eigen::Vector3<T>> c)
{
    return (b - a).cross(c - a);
}

/// @brief Computes the normal vector of a triangle.
/// @param a The first vertex of the triangle.
/// @param b The second vertex of the triangle.
/// @param c The third vertex of the triangle.
/// @return The normal vector of the triangle.
template <typename T>
inline Eigen::Vector3<T> triangle_normal(
    Eigen::ConstRef<Eigen::Vector3<T>> a,
    Eigen::ConstRef<Eigen::Vector3<T>> b,
    Eigen::ConstRef<Eigen::Vector3<T>> c)
{
    return triangle_unnormalized_normal(a, b, c).normalized();
}

/// @brief Computes the Jacobian of the unnormalized normal vector of a
/// triangle.
/// @param a The first vertex of the triangle.
/// @param b The second vertex of the triangle.
/// @param c The third vertex of the triangle.
/// @return The Jacobian of the unnormalized normal vector of the triangle.
template <typename T>
inline Eigen::Matrix<T, 3, 9> triangle_unnormalized_normal_jacobian(
    Eigen::ConstRef<Eigen::Vector3<T>> a,
    Eigen::ConstRef<Eigen::Vector3<T>> b,
    Eigen::ConstRef<Eigen::Vector3<T>> c)
{
    Eigen::Matrix<T, 3, 9> J;
    J.template middleCols<3>(0) = cross_product_matrix<T>(c - b); // ∂n/∂a
    J.template middleCols<3>(3) = cross_product_matrix<T>(a - c); // ∂n/∂b
    J.template middleCols<3>(6) = cross_product_matrix<T>(b - a); // ∂n/∂c
    return J;
}

/// @brief Computes the Hessian of the unnormalized normal vector of a
/// triangle.
/// @param a The first vertex of the triangle.
/// @param b The second vertex of the triangle.
/// @param c The third vertex of the triangle.
/// @return The Hessian of the unnormalized normal vector of the triangle.
template <typename T>
Eigen::Matrix<T, 27, 9> triangle_unnormalized_normal_hessian(
    Eigen::ConstRef<Eigen::Vector3<T>> a,
    Eigen::ConstRef<Eigen::Vector3<T>> b,
    Eigen::ConstRef<Eigen::Vector3<T>> c);

/// @brief Computes the Jacobian of the normal vector of a triangle.
/// @param a The first vertex of the triangle.
/// @param b The second vertex of the triangle.
/// @param c The third vertex of the triangle.
/// @return The Jacobian of the normal vector of the triangle.
template <typename T>
inline Eigen::Matrix<T, 3, 9> triangle_normal_jacobian(
    Eigen::ConstRef<Eigen::Vector3<T>> a,
    Eigen::ConstRef<Eigen::Vector3<T>> b,
    Eigen::ConstRef<Eigen::Vector3<T>> c)
{
    // ∂n̂/∂x = ∂n̂/∂n * ∂n/∂x
    return normalization_jacobian<T>(triangle_unnormalized_normal(a, b, c))
        * triangle_unnormalized_normal_jacobian<T>(a, b, c);
}

/// @brief Computes the Hessian of the normal vector of a triangle.
/// @param a The first vertex of the triangle.
/// @param b The second vertex of the triangle.
/// @param c The third vertex of the triangle.
/// @return The Hessian of the normal vector of the triangle.
template <typename T>
Eigen::Matrix<T, 27, 9> triangle_normal_hessian(
    Eigen::ConstRef<Eigen::Vector3<T>> a,
    Eigen::ConstRef<Eigen::Vector3<T>> b,
    Eigen::ConstRef<Eigen::Vector3<T>> c);

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
template <typename T>
inline Eigen::Vector3<T> line_line_unnormalized_normal(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1)
{
    return (ea1 - ea0).cross(eb1 - eb0);
}

/// @brief Computes the normal vector of two lines.
/// @param ea0 The first vertex of the first line.
/// @param ea1 The second vertex of the first line.
/// @param eb0 The first vertex of the second line.
/// @param eb1 The second vertex of the second line.
/// @return The normal vector of the two lines.
template <typename T>
inline Eigen::Vector3<T> line_line_normal(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1)
{
    return line_line_unnormalized_normal(ea0, ea1, eb0, eb1).normalized();
}

/// @brief Computes the Jacobian of the unnormalized normal vector of two
/// lines.
/// @param ea0 The first vertex of the first line.
/// @param ea1 The second vertex of the first line.
/// @param eb0 The first vertex of the second line.
/// @param eb1 The second vertex of the second line.
/// @return The Jacobian of the unnormalized normal vector of the two lines.
template <typename T>
inline Eigen::Matrix<T, 3, 12> line_line_unnormalized_normal_jacobian(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1)
{
    Eigen::Matrix<T, 3, 12> J;
    J.template middleCols<3>(0) = cross_product_matrix<T>(eb1 - eb0);
    J.template middleCols<3>(3) = cross_product_matrix<T>(eb0 - eb1);
    J.template middleCols<3>(6) = cross_product_matrix<T>(ea0 - ea1);
    J.template middleCols<3>(9) = cross_product_matrix<T>(ea1 - ea0);
    return J;
}

/// @brief Computes the Jacobian of the normal vector of two lines.
/// @param ea0 The first vertex of the first line.
/// @param ea1 The second vertex of the first line.
/// @param eb0 The first vertex of the second line.
/// @param eb1 The second vertex of the second line.
/// @return The Jacobian of the normal vector of the two lines.
template <typename T>
inline Eigen::Matrix<T, 3, 12> line_line_normal_jacobian(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1)
{
    // ∂n̂/∂x = ∂n̂/∂n * ∂n/∂x
    return normalization_jacobian<T>(
               line_line_unnormalized_normal(ea0, ea1, eb0, eb1))
        * line_line_unnormalized_normal_jacobian(ea0, ea1, eb0, eb1);
}

/// @brief Computes the Hessian of the unnormalized normal vector of two
/// lines.
/// @param ea0 The first vertex of the first line.
/// @param ea1 The second vertex of the first line.
/// @param eb0 The first vertex of the second line.
/// @param eb1 The second vertex of the second line.
/// @return The Hessian of the unnormalized normal vector of the two lines.
template <typename T>
Eigen::Matrix<T, 36, 12> line_line_unnormalized_normal_hessian(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1);

/// @brief Computes the Hessian of the normal vector of two lines.
/// @param ea0 The first vertex of the first line.
/// @param ea1 The second vertex of the first line.
/// @param eb0 The first vertex of the second line.
/// @param eb1 The second vertex of the second line.
/// @return The Hessian of the normal vector of the two lines.
template <typename T>
Eigen::Matrix<T, 36, 12> line_line_normal_hessian(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1);

/** @} */

// --- EigenExpression wrappers ---

template <
    typename DerivedX,
    std::enable_if_t<std::is_class_v<DerivedX>, int> = 0>
inline auto normalization_and_jacobian(const Eigen::MatrixBase<DerivedX>& x)
{
    using T = typename DerivedX::Scalar;
    return normalization_and_jacobian(Eigen::Ref<const VectorMax3<T>>(x));
}

template <
    typename DerivedX,
    std::enable_if_t<std::is_class_v<DerivedX>, int> = 0>
inline auto normalization_jacobian(const Eigen::MatrixBase<DerivedX>& x)
{
    using T = typename DerivedX::Scalar;
    return normalization_jacobian(Eigen::Ref<const VectorMax3<T>>(x));
}

template <
    typename DerivedX,
    std::enable_if_t<std::is_class_v<DerivedX>, int> = 0>
inline auto
normalization_and_jacobian_and_hessian(const Eigen::MatrixBase<DerivedX>& x)
{
    using T = typename DerivedX::Scalar;
    return normalization_and_jacobian_and_hessian(
        Eigen::Ref<const VectorMax3<T>>(x));
}

template <
    typename DerivedX,
    std::enable_if_t<std::is_class_v<DerivedX>, int> = 0>
inline auto normalization_hessian(const Eigen::MatrixBase<DerivedX>& x)
{
    using T = typename DerivedX::Scalar;
    return normalization_hessian(Eigen::Ref<const VectorMax3<T>>(x));
}

template <
    typename DerivedX,
    std::enable_if_t<std::is_class_v<DerivedX>, int> = 0>
inline auto cross_product_matrix(const Eigen::MatrixBase<DerivedX>& v)
    -> Eigen::Matrix<typename DerivedX::Scalar, 3, 3>
{
    using T = typename DerivedX::Scalar;
    return cross_product_matrix(Eigen::Ref<const Eigen::Vector3<T>>(v));
}

template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_line_unnormalized_normal(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
    -> VectorMax3<typename DerivedP::Scalar>
{
    using T = typename DerivedP::Scalar;
    return point_line_unnormalized_normal(
        Eigen::Ref<const VectorMax3<T>>(p), Eigen::Ref<const VectorMax3<T>>(e0),
        Eigen::Ref<const VectorMax3<T>>(e1));
}

template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_line_normal(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
    -> VectorMax3<typename DerivedP::Scalar>
{
    using T = typename DerivedP::Scalar;
    return point_line_normal(
        Eigen::Ref<const VectorMax3<T>>(p), Eigen::Ref<const VectorMax3<T>>(e0),
        Eigen::Ref<const VectorMax3<T>>(e1));
}

template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_line_unnormalized_normal_jacobian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
    -> MatrixMax<typename DerivedP::Scalar, 3, 9>
{
    using T = typename DerivedP::Scalar;
    return point_line_unnormalized_normal_jacobian(
        Eigen::Ref<const VectorMax3<T>>(p), Eigen::Ref<const VectorMax3<T>>(e0),
        Eigen::Ref<const VectorMax3<T>>(e1));
}

template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_line_unnormalized_normal_hessian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
    -> MatrixMax<typename DerivedP::Scalar, 27, 9>
{
    using T = typename DerivedP::Scalar;
    return point_line_unnormalized_normal_hessian(
        Eigen::Ref<const VectorMax3<T>>(p), Eigen::Ref<const VectorMax3<T>>(e0),
        Eigen::Ref<const VectorMax3<T>>(e1));
}

template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_line_normal_jacobian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
    -> MatrixMax<typename DerivedP::Scalar, 3, 9>
{
    using T = typename DerivedP::Scalar;
    return point_line_normal_jacobian(
        Eigen::Ref<const VectorMax3<T>>(p), Eigen::Ref<const VectorMax3<T>>(e0),
        Eigen::Ref<const VectorMax3<T>>(e1));
}

template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_line_normal_hessian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
    -> MatrixMax<typename DerivedP::Scalar, 27, 9>
{
    using T = typename DerivedP::Scalar;
    return point_line_normal_hessian(
        Eigen::Ref<const VectorMax3<T>>(p), Eigen::Ref<const VectorMax3<T>>(e0),
        Eigen::Ref<const VectorMax3<T>>(e1));
}

template <typename DerivedA, typename DerivedB, typename DerivedC>
inline auto triangle_unnormalized_normal(
    const Eigen::MatrixBase<DerivedA>& a,
    const Eigen::MatrixBase<DerivedB>& b,
    const Eigen::MatrixBase<DerivedC>& c)
    -> Eigen::Vector3<typename DerivedA::Scalar>
{
    using T = typename DerivedA::Scalar;
    return triangle_unnormalized_normal(
        Eigen::Ref<const Eigen::Vector3<T>>(a),
        Eigen::Ref<const Eigen::Vector3<T>>(b),
        Eigen::Ref<const Eigen::Vector3<T>>(c));
}

template <typename DerivedA, typename DerivedB, typename DerivedC>
inline auto triangle_normal(
    const Eigen::MatrixBase<DerivedA>& a,
    const Eigen::MatrixBase<DerivedB>& b,
    const Eigen::MatrixBase<DerivedC>& c)
    -> Eigen::Vector3<typename DerivedA::Scalar>
{
    using T = typename DerivedA::Scalar;
    return triangle_normal(
        Eigen::Ref<const Eigen::Vector3<T>>(a),
        Eigen::Ref<const Eigen::Vector3<T>>(b),
        Eigen::Ref<const Eigen::Vector3<T>>(c));
}

template <typename DerivedA, typename DerivedB, typename DerivedC>
inline auto triangle_unnormalized_normal_jacobian(
    const Eigen::MatrixBase<DerivedA>& a,
    const Eigen::MatrixBase<DerivedB>& b,
    const Eigen::MatrixBase<DerivedC>& c)
    -> Eigen::Matrix<typename DerivedA::Scalar, 3, 9>
{
    using T = typename DerivedA::Scalar;
    return triangle_unnormalized_normal_jacobian(
        Eigen::Ref<const Eigen::Vector3<T>>(a),
        Eigen::Ref<const Eigen::Vector3<T>>(b),
        Eigen::Ref<const Eigen::Vector3<T>>(c));
}

template <typename DerivedA, typename DerivedB, typename DerivedC>
inline auto triangle_unnormalized_normal_hessian(
    const Eigen::MatrixBase<DerivedA>& a,
    const Eigen::MatrixBase<DerivedB>& b,
    const Eigen::MatrixBase<DerivedC>& c)
    -> Eigen::Matrix<typename DerivedA::Scalar, 27, 9>
{
    using T = typename DerivedA::Scalar;
    return triangle_unnormalized_normal_hessian(
        Eigen::Ref<const Eigen::Vector3<T>>(a),
        Eigen::Ref<const Eigen::Vector3<T>>(b),
        Eigen::Ref<const Eigen::Vector3<T>>(c));
}

template <typename DerivedA, typename DerivedB, typename DerivedC>
inline auto triangle_normal_jacobian(
    const Eigen::MatrixBase<DerivedA>& a,
    const Eigen::MatrixBase<DerivedB>& b,
    const Eigen::MatrixBase<DerivedC>& c)
    -> Eigen::Matrix<typename DerivedA::Scalar, 3, 9>
{
    using T = typename DerivedA::Scalar;
    return triangle_normal_jacobian(
        Eigen::Ref<const Eigen::Vector3<T>>(a),
        Eigen::Ref<const Eigen::Vector3<T>>(b),
        Eigen::Ref<const Eigen::Vector3<T>>(c));
}

template <typename DerivedA, typename DerivedB, typename DerivedC>
inline auto triangle_normal_hessian(
    const Eigen::MatrixBase<DerivedA>& a,
    const Eigen::MatrixBase<DerivedB>& b,
    const Eigen::MatrixBase<DerivedC>& c)
    -> Eigen::Matrix<typename DerivedA::Scalar, 27, 9>
{
    using T = typename DerivedA::Scalar;
    return triangle_normal_hessian(
        Eigen::Ref<const Eigen::Vector3<T>>(a),
        Eigen::Ref<const Eigen::Vector3<T>>(b),
        Eigen::Ref<const Eigen::Vector3<T>>(c));
}

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto line_line_unnormalized_normal(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
    -> Eigen::Vector3<typename DerivedEA0::Scalar>
{
    using T = typename DerivedEA0::Scalar;
    return line_line_unnormalized_normal(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1));
}

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto line_line_normal(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
    -> Eigen::Vector3<typename DerivedEA0::Scalar>
{
    using T = typename DerivedEA0::Scalar;
    return line_line_normal(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1));
}

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto line_line_unnormalized_normal_jacobian(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
    -> Eigen::Matrix<typename DerivedEA0::Scalar, 3, 12>
{
    using T = typename DerivedEA0::Scalar;
    return line_line_unnormalized_normal_jacobian(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1));
}

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto line_line_normal_jacobian(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
    -> Eigen::Matrix<typename DerivedEA0::Scalar, 3, 12>
{
    using T = typename DerivedEA0::Scalar;
    return line_line_normal_jacobian(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1));
}

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto line_line_unnormalized_normal_hessian(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
    -> Eigen::Matrix<typename DerivedEA0::Scalar, 36, 12>
{
    using T = typename DerivedEA0::Scalar;
    return line_line_unnormalized_normal_hessian(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1));
}

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto line_line_normal_hessian(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
    -> Eigen::Matrix<typename DerivedEA0::Scalar, 36, 12>
{
    using T = typename DerivedEA0::Scalar;
    return line_line_normal_hessian(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1));
}

} // namespace ipc
