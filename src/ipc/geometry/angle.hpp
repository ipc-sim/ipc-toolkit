#pragma once

#include <ipc/config.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

namespace detail {
    /// @note Prefer the ipc::dihedral_angle front end.
    template <typename T>
    IPC_TOOLKIT_HOST_DEVICE T dihedral_angle(
        Eigen::ConstRef<Eigen::Vector3<T>> x0,
        Eigen::ConstRef<Eigen::Vector3<T>> x1,
        Eigen::ConstRef<Eigen::Vector3<T>> x2,
        Eigen::ConstRef<Eigen::Vector3<T>> x3);

    /// @note Prefer the ipc::dihedral_angle_gradient front end.
    template <typename T>
    IPC_TOOLKIT_HOST_DEVICE Eigen::Vector<T, 12> dihedral_angle_gradient(
        Eigen::ConstRef<Eigen::Vector3<T>> x0,
        Eigen::ConstRef<Eigen::Vector3<T>> x1,
        Eigen::ConstRef<Eigen::Vector3<T>> x2,
        Eigen::ConstRef<Eigen::Vector3<T>> x3);

    /// @note Prefer the ipc::dihedral_angle_hessian front end.
    template <typename T>
    IPC_TOOLKIT_HOST_DEVICE Eigen::Matrix<T, 12, 12> dihedral_angle_hessian(
        Eigen::ConstRef<Eigen::Vector3<T>> x0,
        Eigen::ConstRef<Eigen::Vector3<T>> x1,
        Eigen::ConstRef<Eigen::Vector3<T>> x2,
        Eigen::ConstRef<Eigen::Vector3<T>> x3);
} // namespace detail

/// @brief Compute the bending angle between two triangles sharing an edge.
///     x0---x3
///      | \ |
///     x2---x1
/// @param x0 The first vertex of the edge.
/// @param x1 The second vertex of the edge.
/// @param x2 The opposite vertex of the first triangle.
/// @param x3 The opposite vertex of the second triangle.
/// @return The bending angle between the two triangles.
template <
    typename DerivedX0,
    typename DerivedX1,
    typename DerivedX2,
    typename DerivedX3>
IPC_TOOLKIT_HOST_DEVICE inline auto dihedral_angle(
    const Eigen::MatrixBase<DerivedX0>& x0,
    const Eigen::MatrixBase<DerivedX1>& x1,
    const Eigen::MatrixBase<DerivedX2>& x2,
    const Eigen::MatrixBase<DerivedX3>& x3)
{
    using T = typename DerivedX0::Scalar;
    return detail::dihedral_angle<T>(x0, x1, x2, x3);
}

/// @brief Compute the Jacobian of the bending angle between two triangles sharing an edge.
///     x0---x3
///      | \ |
///     x2---x1
/// @param x0 The first vertex of the edge.
/// @param x1 The second vertex of the edge.
/// @param x2 The opposite vertex of the first triangle.
/// @param x3 The opposite vertex of the second triangle.
/// @return The Jacobian matrix of the bending angle with respect to the input vertices.
template <
    typename DerivedX0,
    typename DerivedX1,
    typename DerivedX2,
    typename DerivedX3>
IPC_TOOLKIT_HOST_DEVICE inline auto dihedral_angle_gradient(
    const Eigen::MatrixBase<DerivedX0>& x0,
    const Eigen::MatrixBase<DerivedX1>& x1,
    const Eigen::MatrixBase<DerivedX2>& x2,
    const Eigen::MatrixBase<DerivedX3>& x3)
{
    using T = typename DerivedX0::Scalar;
    return detail::dihedral_angle_gradient<T>(x0, x1, x2, x3);
}

/// @brief Compute the Hessian of the bending angle between two triangles sharing an edge.
///     x0---x3
///      | \ |
///     x2---x1
/// @param x0 The first vertex of the edge.
/// @param x1 The second vertex of the edge.
/// @param x2 The opposite vertex of the first triangle.
/// @param x3 The opposite vertex of the second triangle.
/// @return The 12x12 Hessian matrix of the bending angle with respect to the input vertices.
template <
    typename DerivedX0,
    typename DerivedX1,
    typename DerivedX2,
    typename DerivedX3>
IPC_TOOLKIT_HOST_DEVICE inline auto dihedral_angle_hessian(
    const Eigen::MatrixBase<DerivedX0>& x0,
    const Eigen::MatrixBase<DerivedX1>& x1,
    const Eigen::MatrixBase<DerivedX2>& x2,
    const Eigen::MatrixBase<DerivedX3>& x3)
{
    using T = typename DerivedX0::Scalar;
    return detail::dihedral_angle_hessian<T>(x0, x1, x2, x3);
}

} // namespace ipc
