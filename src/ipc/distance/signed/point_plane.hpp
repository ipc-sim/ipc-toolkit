#pragma once

#include <ipc/geometry/normal.hpp>

namespace ipc {

/// Compute the signed distance from a point to the plane of a triangle.
///
/// The signed distance is computed as:
///   d = triangle_normal(t0, t1, t2) · (p - t0)
/// The sign corresponds to the orientation of the triangle normal.
///
/// @param p   The query point (3D).
/// @param t0  First vertex of the triangle (3D), used as a point on the plane.
/// @param t1  Second vertex of the triangle (3D).
/// @param t2  Third vertex of the triangle (3D).
/// @return    The signed distance from p to the plane of the triangle.
inline double point_plane_signed_distance(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> t0,
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2)
{
    return triangle_normal(t0, t1, t2).dot(p - t0);
}

/// Compute the gradient of the signed point-to-plane distance.
///
/// The returned Vector12d contains derivatives in the order:
///   [ ∂d/∂p (3), ∂d/∂t0 (3), ∂d/∂t1 (3), ∂d/∂t2 (3) ]
///
/// @param p   The query point (3D).
/// @param t0  First vertex of the triangle (3D).
/// @param t1  Second vertex of the triangle (3D).
/// @param t2  Third vertex of the triangle (3D).
/// @return    A Vector12d containing the gradient of the signed distance.
Vector12d point_plane_signed_distance_gradient(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> t0,
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2);

/// Compute the Hessian (matrix of second derivatives) of the signed distance.
///
/// The returned Matrix12d is the 12x12 Hessian with variables ordered as:
///   [ p (3), t0 (3), t1 (3), t2 (3) ].
///
/// @param p   The query point (3D).
/// @param t0  First vertex of the triangle (3D).
/// @param t1  Second vertex of the triangle (3D).
/// @param t2  Third vertex of the triangle (3D).
/// @return    A Matrix12d representing the Hessian of the signed distance.
Matrix12d point_plane_signed_distance_hessian(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> t0,
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2);

} // namespace ipc