#pragma once

#include <ipc/geometry/normal.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

namespace detail {
    /// Compute the signed distance between two lines in 3D.
    ///
    /// The two lines are specified by points (ea0, ea1) and (eb0, eb1).
    /// The returned value is the scalar projection of the vector (ea0 - eb0)
    /// onto the unit normal that is perpendicular to both lines (as produced by
    /// line_line_normal). In other words, this is the signed distance along the
    /// common normal between the two lines; the sign follows the orientation of
    /// the normal returned by line_line_normal.
    ///
    /// @param ea0 First endpoint (or a point) on the first line.
    /// @param ea1 Second endpoint (or a direction-defining point) on the first line.
    /// @param eb0 First endpoint (or a point) on the second line.
    /// @param eb1 Second endpoint (or a direction-defining point) on the second line.
    /// @return Signed distance between the two lines along the common normal. If
    /// the lines are parallel (no unique common normal), the result may be
    /// ill-conditioned or implementation-defined.
    /// @note The points may be any two distinct points on each line (they need not
    /// represent segment endpoints). Behavior is undefined if ea0 == ea1 or eb0
    /// == eb1.
    template <typename T>
    inline T line_line_signed_distance(
        Eigen::ConstRef<Eigen::Vector3<T>> ea0,
        Eigen::ConstRef<Eigen::Vector3<T>> ea1,
        Eigen::ConstRef<Eigen::Vector3<T>> eb0,
        Eigen::ConstRef<Eigen::Vector3<T>> eb1)
    {
        return line_line_normal(ea0, ea1, eb0, eb1).dot(ea0 - eb0);
    }

    /// Compute the gradient of the signed line-line distance with respect to
    /// the four 3D input points: ea0, ea1, eb0, eb1.
    ///
    /// The returned gradient is a 12-vector ordered as [d/d(ea0); d/d(ea1);
    /// d/d(eb0); d/d(eb1)], where each block is a 3-vector. This function
    /// differentiates the signed distance produced by
    /// line_line_signed_distance. The gradient may be undefined or numerically
    /// unstable when the two lines are parallel or when the direction-defining
    /// points coincide.
    ///
    /// @param ea0 First point on the first line (as in line_line_signed_distance).
    /// @param ea1 Second point on the first line.
    /// @param eb0 First point on the second line.
    /// @param eb1 Second point on the second line.
    /// @return 12-dimensional gradient of the signed distance with respect to
    /// [ea0, ea1, eb0, eb1].
    /// @pre ea0 != ea1 and eb0 != eb1. For near-parallel lines, results may be unstable.
    ///
    /// @see line_line_signed_distance, line_line_normal
    template <typename T>
    inline Eigen::Vector<T, 12> line_line_signed_distance_gradient(
        Eigen::ConstRef<Eigen::Vector3<T>> ea0,
        Eigen::ConstRef<Eigen::Vector3<T>> ea1,
        Eigen::ConstRef<Eigen::Vector3<T>> eb0,
        Eigen::ConstRef<Eigen::Vector3<T>> eb1)
    {
        const Eigen::Vector3<T> n = line_line_normal(ea0, ea1, eb0, eb1);
        const Eigen::Matrix<T, 3, 12> jac_n =
            line_line_normal_jacobian(ea0, ea1, eb0, eb1);
        Eigen::Vector<T, 12> grad = jac_n.transpose() * (ea0 - eb0);
        grad.template segment<3>(0) += n;
        grad.template segment<3>(6) -= n;
        return grad;
    }

    /// Compute the Hessian (second derivative) of the signed line-line distance
    /// with respect to the four 3D input points: ea0, ea1, eb0, eb1.
    ///
    /// The returned matrix is 12x12 and corresponds to the second derivatives
    /// of the scalar signed distance with respect to the stacked variable
    /// [ea0; ea1; eb0; eb1]. The Hessian captures curvature information and is
    /// typically required for second-order optimization or sensitivity
    /// analysis. As with the gradient, the Hessian is undefined or numerically
    /// unstable if either input line is degenerate (coincident points) or if
    /// the lines are parallel (no unique perpendicular direction).
    ///
    /// @param ea0 First point on the first line.
    /// @param ea1 Second point on the first line.
    /// @param eb0 First point on the second line.
    /// @param eb1 Second point on the second line.
    /// @return 12x12 Hessian matrix of second derivatives of the signed distance.
    /// @pre ea0 != ea1 and eb0 != eb1. Results may be invalid for parallel lines.
    ///
    /// @see line_line_signed_distance, line_line_signed_distance_gradient
    template <typename T>
    Eigen::Matrix<T, 12, 12> line_line_signed_distance_hessian(
        Eigen::ConstRef<Eigen::Vector3<T>> ea0,
        Eigen::ConstRef<Eigen::Vector3<T>> ea1,
        Eigen::ConstRef<Eigen::Vector3<T>> eb0,
        Eigen::ConstRef<Eigen::Vector3<T>> eb1);
} // namespace detail

/// Compute the signed distance between two lines in 3D.
///
/// The two lines are specified by points (ea0, ea1) and (eb0, eb1).
/// The returned value is the scalar projection of the vector (ea0 - eb0)
/// onto the unit normal that is perpendicular to both lines (as produced by
/// line_line_normal). In other words, this is the signed distance along the
/// common normal between the two lines; the sign follows the orientation of
/// the normal returned by line_line_normal.
///
/// @param ea0 First endpoint (or a point) on the first line.
/// @param ea1 Second endpoint (or a direction-defining point) on the first line.
/// @param eb0 First endpoint (or a point) on the second line.
/// @param eb1 Second endpoint (or a direction-defining point) on the second line.
/// @return Signed distance between the two lines along the common normal. If
/// the lines are parallel (no unique common normal), the result may be
/// ill-conditioned or implementation-defined.
/// @note The points may be any two distinct points on each line (they need not
/// represent segment endpoints). Behavior is undefined if ea0 == ea1 or eb0
/// == eb1.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto line_line_signed_distance(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
{
    using T = typename DerivedEA0::Scalar;
    // NOTE: explicit <T>: the detail overload cannot deduce T from an
    // arbitrary Eigen expression.
    return detail::line_line_signed_distance<T>(ea0, ea1, eb0, eb1);
}

/// Compute the gradient of the signed line-line distance with respect to
/// the four 3D input points: ea0, ea1, eb0, eb1.
///
/// The returned gradient is a 12-vector ordered as [d/d(ea0); d/d(ea1);
/// d/d(eb0); d/d(eb1)], where each block is a 3-vector. This function
/// differentiates the signed distance produced by
/// line_line_signed_distance. The gradient may be undefined or numerically
/// unstable when the two lines are parallel or when the direction-defining
/// points coincide.
///
/// @param ea0 First point on the first line (as in line_line_signed_distance).
/// @param ea1 Second point on the first line.
/// @param eb0 First point on the second line.
/// @param eb1 Second point on the second line.
/// @return 12-dimensional gradient of the signed distance with respect to
/// [ea0, ea1, eb0, eb1].
/// @pre ea0 != ea1 and eb0 != eb1. For near-parallel lines, results may be unstable.
///
/// @see line_line_signed_distance, line_line_normal
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto line_line_signed_distance_gradient(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
{
    using T = typename DerivedEA0::Scalar;
    return detail::line_line_signed_distance_gradient<T>(ea0, ea1, eb0, eb1);
}

/// Compute the Hessian (second derivative) of the signed line-line distance
/// with respect to the four 3D input points: ea0, ea1, eb0, eb1.
///
/// The returned matrix is 12x12 and corresponds to the second derivatives
/// of the scalar signed distance with respect to the stacked variable
/// [ea0; ea1; eb0; eb1]. The Hessian captures curvature information and is
/// typically required for second-order optimization or sensitivity
/// analysis. As with the gradient, the Hessian is undefined or numerically
/// unstable if either input line is degenerate (coincident points) or if
/// the lines are parallel (no unique perpendicular direction).
///
/// @param ea0 First point on the first line.
/// @param ea1 Second point on the first line.
/// @param eb0 First point on the second line.
/// @param eb1 Second point on the second line.
/// @return 12x12 Hessian matrix of second derivatives of the signed distance.
/// @pre ea0 != ea1 and eb0 != eb1. Results may be invalid for parallel lines.
///
/// @see line_line_signed_distance, line_line_signed_distance_gradient
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto line_line_signed_distance_hessian(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
{
    using T = typename DerivedEA0::Scalar;
    return detail::line_line_signed_distance_hessian<T>(ea0, ea1, eb0, eb1);
}

} // namespace ipc
