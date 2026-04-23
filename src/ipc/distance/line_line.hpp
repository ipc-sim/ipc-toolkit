#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Geometry>

namespace ipc {

// Symbolically generated derivatives
namespace autogen {
    // clang-format off
    template <typename T>
    void line_line_distance_gradient(
        T v01, T v02, T v03, T v11, T v12, T v13, T v21, T v22, T v23, T v31, T v32, T v33, T g[12]);
    template <typename T>
    void line_line_distance_hessian(
        T v01, T v02, T v03, T v11, T v12, T v13, T v21, T v22, T v23, T v31, T v32, T v33, T H[144]);
    // clang-format on
} // namespace autogen

/// @brief Compute the distance between a two infinite lines in 3D.
/// @note The distance is actually squared distance.
/// @warning If the lines are parallel this function returns a distance of zero.
/// @param ea0 The first vertex of the edge defining the first line.
/// @param ea1 The second vertex of the edge defining the first line.
/// @param eb0 The first vertex of the edge defining the second line.
/// @param eb1 The second vertex of the edge defining the second line.
/// @return The distance between the two lines.
template <typename T>
inline T line_line_distance(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1)
{
    const Eigen::Vector3<T> normal = (ea1 - ea0).cross(eb1 - eb0);
    const T line_to_line = (eb0 - ea0).dot(normal);
    return line_to_line * line_to_line / normal.squaredNorm();
}

/// @brief Compute the gradient of the distance between a two lines in 3D.
/// @note The distance is actually squared distance.
/// @warning If the lines are parallel this function returns a distance of zero.
/// @param ea0 The first vertex of the edge defining the first line.
/// @param ea1 The second vertex of the edge defining the first line.
/// @param eb0 The first vertex of the edge defining the second line.
/// @param eb1 The second vertex of the edge defining the second line.
/// @return The gradient of the distance wrt ea0, ea1, eb0, and eb1.
template <typename T>
inline Eigen::Vector<T, 12> line_line_distance_gradient(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1)
{
    Eigen::Vector<T, 12> grad;
    autogen::line_line_distance_gradient(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], grad.data());
    return grad;
}

/// @brief Compute the hessian of the distance between a two lines in 3D.
/// @note The distance is actually squared distance.
/// @warning If the lines are parallel this function returns a distance of zero.
/// @param ea0 The first vertex of the edge defining the first line.
/// @param ea1 The second vertex of the edge defining the first line.
/// @param eb0 The first vertex of the edge defining the second line.
/// @param eb1 The second vertex of the edge defining the second line.
/// @return The hessian of the distance wrt ea0, ea1, eb0, and eb1.
template <typename T>
inline Eigen::Matrix<T, 12, 12> line_line_distance_hessian(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1)
{
    Eigen::Matrix<T, 12, 12> hess;
    autogen::line_line_distance_hessian(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], hess.data());
    return hess;
}

// --- EigenExpression wrappers ---

template <
    EigenExpression DerivedEa0,
    EigenExpression DerivedEa1,
    EigenExpression DerivedEb0,
    EigenExpression DerivedEb1>
inline auto line_line_distance(
    const Eigen::MatrixBase<DerivedEa0>& ea0,
    const Eigen::MatrixBase<DerivedEa1>& ea1,
    const Eigen::MatrixBase<DerivedEb0>& eb0,
    const Eigen::MatrixBase<DerivedEb1>& eb1) -> typename DerivedEa0::Scalar
{
    using T = typename DerivedEa0::Scalar;
    return line_line_distance(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1));
}

template <
    EigenExpression DerivedEa0,
    EigenExpression DerivedEa1,
    EigenExpression DerivedEb0,
    EigenExpression DerivedEb1>
inline auto line_line_distance_gradient(
    const Eigen::MatrixBase<DerivedEa0>& ea0,
    const Eigen::MatrixBase<DerivedEa1>& ea1,
    const Eigen::MatrixBase<DerivedEb0>& eb0,
    const Eigen::MatrixBase<DerivedEb1>& eb1)
    -> Eigen::Vector<typename DerivedEa0::Scalar, 12>
{
    using T = typename DerivedEa0::Scalar;
    return line_line_distance_gradient(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1));
}

template <
    EigenExpression DerivedEa0,
    EigenExpression DerivedEa1,
    EigenExpression DerivedEb0,
    EigenExpression DerivedEb1>
inline auto line_line_distance_hessian(
    const Eigen::MatrixBase<DerivedEa0>& ea0,
    const Eigen::MatrixBase<DerivedEa1>& ea1,
    const Eigen::MatrixBase<DerivedEb0>& eb0,
    const Eigen::MatrixBase<DerivedEb1>& eb1)
    -> Eigen::Matrix<typename DerivedEa0::Scalar, 12, 12>
{
    using T = typename DerivedEa0::Scalar;
    return line_line_distance_hessian(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1));
}

} // namespace ipc
