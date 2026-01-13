#pragma once

#include <ipc/geometry/normal.hpp>

namespace ipc {

/// Compute the signed distance from a point to a directed line segment.
///
/// The signed distance is computed as d = n · (p - e0),
/// where n = point_line_normal(p, e0, e1) is the unit normal associated with
/// the directed edge (e0 -> e1) chosen consistently for the point p.
/// Positive/negative sign indicates which side of the directed edge the point
/// lies on.
///
/// @param p   The query point (2D).
/// @param e0  The first endpoint of the directed edge (2D).
/// @param e1  The second endpoint of the directed edge (2D).
/// @return    The signed scalar distance from p to the infinite line through e0 and e1.
/// @note      The edge must be non-degenerate (e0 != e1).
inline double point_line_signed_distance(
    Eigen::ConstRef<Eigen::Vector2d> p,
    Eigen::ConstRef<Eigen::Vector2d> e0,
    Eigen::ConstRef<Eigen::Vector2d> e1)
{
    return point_line_normal(p, e0, e1).dot(p - e0);
}

/// Compute the gradient of the signed point-to-line distance with respect to
/// all input coordinates packed as [p_x, p_y, e0_x, e0_y, e1_x, e1_y]^T.
///
/// The returned Vector6d contains partial derivatives of the scalar signed
/// distance d with respect to each coordinate in the above ordering:
///   grad = [∂d/∂p, ∂d/∂e0, ∂d/∂e1]^T.
///
/// @param p   The query point (2D).
/// @param e0  The first endpoint of the directed edge (2D).
/// @param e1  The second endpoint of the directed edge (2D).
/// @return    A 6-vector containing the gradient of the signed distance.
/// @note      The edge must be non-degenerate (e0 != e1).
Vector6d point_line_signed_distance_gradient(
    Eigen::ConstRef<Eigen::Vector2d> p,
    Eigen::ConstRef<Eigen::Vector2d> e0,
    Eigen::ConstRef<Eigen::Vector2d> e1);

/// Compute the Hessian (second derivatives) of the signed point-to-line
/// distance with respect to all input coordinates packed as [p_x, p_y, e0_x,
/// e0_y, e1_x, e1_y]^T.
///
/// The returned Matrix6d is the symmetric 6x6 matrix of second partial
/// derivatives:
///   H_ij = ∂^2 d / (∂x_i ∂x_j),
/// with the same coordinate ordering as in the gradient.
///
/// @param p   The query point (2D).
/// @param e0  The first endpoint of the directed edge (2D).
/// @param e1  The second endpoint of the directed edge (2D).
/// @return    A 6x6 Hessian matrix of the signed distance.
/// @note      The edge must be non-degenerate (e0 != e1).
Matrix6d point_line_signed_distance_hessian(
    Eigen::ConstRef<Eigen::Vector2d> p,
    Eigen::ConstRef<Eigen::Vector2d> e0,
    Eigen::ConstRef<Eigen::Vector2d> e1);

} // namespace ipc