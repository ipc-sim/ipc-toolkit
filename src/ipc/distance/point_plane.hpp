#pragma once

#include <ipc/config.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Compute the distance between a point and a plane.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param origin The origin of the plane.
/// @param normal The normal of the plane.
/// @return The distance between the point and plane.
IPC_TOOLKIT_INLINE double point_plane_distance(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> origin,
    Eigen::ConstRef<Eigen::Vector3d> normal);

/// @brief Compute the distance between a point and a plane.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param t0 The first vertex of the triangle.
/// @param t1 The second vertex of the triangle.
/// @param t2 The third vertex of the triangle.
/// @return The distance between the point and plane.
IPC_TOOLKIT_INLINE double point_plane_distance(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> t0,
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2);

/// @brief Compute the gradient of the distance between a point and a plane.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param origin The origin of the plane.
/// @param normal The normal of the plane.
/// @return The gradient of the distance wrt p.
IPC_TOOLKIT_INLINE Eigen::Vector3d point_plane_distance_gradient(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> origin,
    Eigen::ConstRef<Eigen::Vector3d> normal);

/// @brief Compute the gradient of the distance between a point and a plane.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param t0 The first vertex of the triangle.
/// @param t1 The second vertex of the triangle.
/// @param t2 The third vertex of the triangle.
/// @return The gradient of the distance wrt p, t0, t1, and t2.
IPC_TOOLKIT_INLINE Vector12d point_plane_distance_gradient(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> t0,
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2);

/// @brief Compute the hessian of the distance between a point and a plane.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param origin The origin of the plane.
/// @param normal The normal of the plane.
/// @return The hessian of the distance wrt p.
IPC_TOOLKIT_INLINE Eigen::Matrix3d point_plane_distance_hessian(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> origin,
    Eigen::ConstRef<Eigen::Vector3d> normal);

/// @brief Compute the hessian of the distance between a point and a plane.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param t0 The first vertex of the triangle.
/// @param t1 The second vertex of the triangle.
/// @param t2 The third vertex of the triangle.
/// @return The hessian of the distance wrt p, t0, t1, and t2.
IPC_TOOLKIT_INLINE Matrix12d point_plane_distance_hessian(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> t0,
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2);

// Symbolically generated derivatives;
namespace autogen {
    IPC_TOOLKIT_INLINE void point_plane_distance_gradient(
        double v01,
        double v02,
        double v03,
        double v11,
        double v12,
        double v13,
        double v21,
        double v22,
        double v23,
        double v31,
        double v32,
        double v33,
        double g[12]);

    IPC_TOOLKIT_INLINE void point_plane_distance_hessian(
        double v01,
        double v02,
        double v03,
        double v11,
        double v12,
        double v13,
        double v21,
        double v22,
        double v23,
        double v31,
        double v32,
        double v33,
        double H[144]);
} // namespace autogen

} // namespace ipc

#ifdef IPC_TOOLKIT_WITH_CUDA
#include "point_plane.cpp"
#endif
