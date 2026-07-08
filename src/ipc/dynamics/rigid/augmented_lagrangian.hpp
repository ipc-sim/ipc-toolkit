#pragma once

#include <ipc/dynamics/rigid/rigid_potential.hpp>

namespace ipc::rigid {

/// @brief Augmented Lagrangian driving KINEMATIC bodies to per-step target
/// poses [Ferguson et al. 2021].
///
/// The equality constraint x = x̂ (per kinematic body) is enforced softly so
/// the barrier can keep the trajectory intersection-free while the body
/// approaches its target. Per kinematic body, with penalties κ_q/κ_Q and
/// multipliers λ ∈ R³ and Λ:
/// \f[
///     E = \tfrac{\kappa_q}{2} m \|q - \hat{q}\|^2
///       - \sqrt{m}\, \lambda^\top (q - \hat{q})
///       + \tfrac{\kappa_Q}{2} \operatorname{tr}((Q - \hat{Q}) J (Q -
///       \hat{Q})^\top) - \operatorname{tr}(\Lambda^\top (Q - \hat{Q})
///       J^{1/2})
/// \f]
/// The mass/inertia weighting makes the penalties scale-invariant. The
/// classic method-of-multipliers schedule (per linear/angular channel):
/// satisfied when the progress η ≥ satisfied_progress (freeze the DOFs);
/// double the penalty while η < stall_progress and κ < max_penalty;
/// otherwise update the multipliers.
///
/// This is a pure potential: the simulator applies the (βΔt)² acceleration
/// scaling at summation (a positive constant scale that does not change the
/// constrained minimizer).
class AugmentedLagrangian : public RigidPotential {
public:
    struct Params {
        /// @brief Initial penalty κ_q = κ_Q at the start of each step.
        double initial_penalty = 1e3;
        /// @brief Maximum penalty (stop doubling).
        double max_penalty = 1e8;
        /// @brief Progress η at which a channel is satisfied (DOFs freeze).
        double satisfied_progress = 0.999;
        /// @brief Progress below which the penalty is doubled (otherwise the
        /// multipliers are updated).
        double stall_progress = 0.99;
    };

    AugmentedLagrangian() = default;

    explicit AugmentedLagrangian(const Params& params) : m_params(params) { }

    /// @brief Start a new step: capture the start poses and targets, reset
    /// the penalties, multipliers, and satisfied flags.
    /// @param bodies The collection of rigid bodies.
    /// @param x0 Step-start DOFs (ndof per body).
    /// @param targets Per-body target poses (only KINEMATIC entries used).
    void init(
        const RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x0,
        const std::vector<Pose>& targets);

    /// @brief Whether any kinematic channel is still unsatisfied.
    bool active() const
    {
        return m_num_kinematic > 0
            && !(m_linear_satisfied && m_angular_satisfied);
    }

    /// @brief Linear progress η_q ∈ (-∞, 1] toward the targets.
    double
    linear_progress(const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x) const;

    /// @brief Angular progress η_Q ∈ (-∞, 1] toward the targets.
    double angular_progress(
        const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x) const;

    /// @brief Apply the update policy at the (energy-converged) state x.
    void update(const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x);

    bool linear_satisfied() const { return m_linear_satisfied; }
    bool angular_satisfied() const { return m_angular_satisfied; }

    double linear_penalty() const { return m_kappa_q; }
    double angular_penalty() const { return m_kappa_Q; }

    const std::vector<VectorMax3d>& linear_multipliers() const
    {
        return m_lambda;
    }
    const std::vector<MatrixMax3d>& angular_multipliers() const
    {
        return m_Lambda;
    }

    // ---- Cumulative functions -----------------------------------------------

    using RigidPotential::operator();
    using RigidPotential::gradient;
    using RigidPotential::hessian;

    // ---- Per-body functions -------------------------------------------------

    /// @brief AL energy of one body (zero unless KINEMATIC and active()).
    double operator()(
        const size_t body_id,
        const RigidBody& body,
        Eigen::ConstRef<VectorMax6d> x) const override;

    VectorMax6d gradient(
        const size_t body_id,
        const RigidBody& body,
        Eigen::ConstRef<VectorMax6d> x) const override;

    MatrixMax6d hessian(
        const size_t body_id,
        const RigidBody& body,
        Eigen::ConstRef<VectorMax6d> x,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const override;

private:
    Params m_params;

    /// @brief Number of KINEMATIC bodies at init.
    size_t m_num_kinematic = 0;

    /// @brief Step-start poses (η denominators).
    std::vector<Pose> m_start;
    /// @brief Target poses.
    std::vector<Pose> m_targets;
    /// @brief Target rotation matrices (cached from m_targets).
    std::vector<MatrixMax3d> m_target_rotations;

    /// @brief Penalties (κ_q linear, κ_Q angular).
    double m_kappa_q = 1e3;
    double m_kappa_Q = 1e3; // NOLINT(readability-identifier-naming)

    /// @brief Per-body multipliers (λ linear, Λ angular).
    std::vector<VectorMax3d> m_lambda;
    std::vector<MatrixMax3d> m_Lambda; // NOLINT(readability-identifier-naming)

    /// @brief Per-channel satisfied flags (freeze the channel's DOFs).
    bool m_linear_satisfied = false;
    bool m_angular_satisfied = false;
};

} // namespace ipc::rigid
