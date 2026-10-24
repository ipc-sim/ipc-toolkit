#pragma once

#include <ipc/ccd/nonlinear_ccd.hpp>
#include <ipc/dynamics/rigid/pose.hpp>

#include <algorithm>

namespace ipc::rigid {

/// @brief Bound the deviation of a rigidly-carried point from the chord of its
/// trajectory.
///
/// Consider a point at body-frame position x̄ carried by a rigid motion with
/// rotation vector interpolated linearly: f(s) = R(θa + s δθ) x̄, s ∈ [0, 1].
/// The deviation from the endpoint linear interpolant l(s) = (1-s) f(0) + s
/// f(1) is bounded by the classic C² interpolation-error bound
///     ‖f(s) − l(s)‖ ≤ ⅛ max_s ‖f″(s)‖.
/// Using the Fréchet-integral representation of d/ds exp([·]×) — where each
/// factor is an orthogonal matrix (spectral norm 1) times [δθ]× (spectral norm
/// ‖δθ‖) — gives ‖f″(s)‖ ≤ ‖δθ‖² ‖x̄‖ *without any fixed-axis assumption*, so
/// the bound is conservative even when θa and θa + δθ are not parallel. In 2D
/// the axis is fixed and this is the familiar sagitta bound
/// r (1 − cos(δ/2)) ≤ r δ²/8. Since both the curve and its chord lie in the
/// ball of radius ‖x̄‖ about the rotation center, the deviation never exceeds
/// 2‖x̄‖, which caps the bound for large rotations.
///
/// See also: Redon et al. [2002], "Fast Continuous Collision Detection between
/// Rigid Bodies."
///
/// @param radius Distance ‖x̄‖ of the point from the rotation center
/// @param dtheta Norm of the rotation-vector change ‖δθ‖ over the interval
/// @return Conservative bound on max_s ‖f(s) − l(s)‖
inline double rigid_chord_deviation(const double radius, const double dtheta)
{
    return std::min(radius * dtheta * dtheta / 8, 2 * radius);
}

/// @brief Bound the chord deviation of a point on body j expressed in body i's
/// rest frame.
///
/// The relative curve is x(t) = Rᵢ(t)ᵀ (Rⱼ(t) x̄ + Δq(t)) with Δq(t) = pⱼ(t) −
/// pᵢ(t) affine in t. Applying the same ⅛ max‖x″‖ interpolation-error bound as
/// rigid_chord_deviation with the product rule x″ = A″u + 2A′u′ + Au″ (A =
/// Rᵢᵀ, u = Rⱼ x̄ + Δq) and the Fréchet bounds ‖A′‖ ≤ ‖δθᵢ‖, ‖A″‖ ≤ ‖δθᵢ‖²
/// gives
///     dev = ⅛ [ (‖δθᵢ‖ + ‖δθⱼ‖)² ‖x̄‖
///               + ‖δθᵢ‖ (‖δθᵢ‖ max‖Δq‖ + 2 ‖Δq(1) − Δq(0)‖) ].
/// Both the curve and its chord lie in the ball of radius ‖x̄‖ + max‖Δq‖ about
/// the origin, capping the deviation at twice that. When body i is static
/// (δθᵢ = 0) this degenerates to the world-space rigid_chord_deviation bound.
///
/// @param radius Distance ‖x̄‖ of the point from body j's rotation center
/// @param dtheta_i Norm of body i's rotation-vector change over the interval
/// @param dtheta_j Norm of body j's rotation-vector change over the interval
/// @param max_translation Maximum of ‖Δq‖ over the interval (attained at an endpoint since Δq is affine)
/// @param translation_change Norm of the relative translation change ‖Δq(1) − Δq(0)‖ over the interval
/// @return Conservative bound on the deviation of the relative curve from its chord
inline double relative_rigid_chord_deviation(
    const double radius,
    const double dtheta_i,
    const double dtheta_j,
    const double max_translation,
    const double translation_change)
{
    const double dtheta_sum = dtheta_i + dtheta_j;
    const double dev =
        (dtheta_sum * dtheta_sum * radius
         + dtheta_i * (dtheta_i * max_translation + 2 * translation_change))
        / 8;
    return std::min(dev, 2 * (radius + max_translation));
}

/// @brief The trajectory of a point carried by a rigid body.
///
/// The body moves from pose_t0 to pose_t1 with position and rotation vector
/// interpolated linearly (matching the incremental potential's pose
/// parameterization):
///     x(t) = R((1−t) θ₀ + t θ₁) x̄ + (1−t) p₀ + t p₁.
class RigidTrajectory : virtual public NonlinearTrajectory {
public:
    /// @brief Construct the trajectory of a point on a rigid body.
    /// @param rest_position The body-frame (rest) position x̄ of the point
    /// @param pose_t0 The pose of the body at the start of the time step
    /// @param pose_t1 The pose of the body at the end of the time step
    RigidTrajectory(
        Eigen::ConstRef<VectorMax3d> rest_position,
        const Pose& pose_t0,
        const Pose& pose_t1);

    /// @brief Compute the point's position at time t ∈ [0, 1].
    VectorMax3d operator()(const double t) const override;

    /// @brief Conservative bound on the deviation from the linearized trajectory.
    ///
    /// The translational part is affine in t, so it matches its chord exactly
    /// and only the rotational part contributes (see rigid_chord_deviation).
    ///
    /// @param t0 Start time of the sub-interval
    /// @param t1 End time of the sub-interval
    double
    max_distance_from_linear(const double t0, const double t1) const override;

private:
    VectorMax3d m_rest_position;
    Pose m_pose_t0;
    Pose m_pose_t1;
    double m_radius;      ///< ‖x̄‖ (cached)
    double m_dtheta_norm; ///< ‖θ₁ − θ₀‖ (cached)
};

} // namespace ipc::rigid
