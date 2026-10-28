#pragma once

#include <ipc/dynamics/rigid/rigid_bodies.hpp>

#include <limits>

namespace ipc::affine {

/// @brief Augmented Lagrangian driving KINEMATIC bodies to per-step target
/// poses in affine coordinates [Ferguson et al. 2021]. The rigid
/// parameterization reaches it through the to-affine map's chain rule.
///
/// With y = [p; vec(A)] per body and the exact ABD mass blocks
/// M = blkdiag(m I₃, W) with W = J ⊗ I₃ (diagonal), the two channels are
/// \f[
///     E = \tfrac{\kappa_q}{2} m \|q - \hat{q}\|^2
///       - s\, \lambda^\top (q - \hat{q})
///       + \tfrac{\kappa_Q}{2} (a - \hat{a})^\top W (a - \hat{a})
///       - s\, \lambda_A^\top (a - \hat{a}),
/// \f]
/// where â = vec(R(θ̂)) embeds the rigid target rotation and s = (βΔt)² is
/// the time integrator's acceleration scaling. Fully quadratic with constant
/// PSD Hessian. We use Rigid IPC's method-of-multipliers schedule, and also
/// double the penalty when progress stagnates short of the target.
///
/// The two halves are scaled differently on purpose. The penalty is a
/// numerical device, so we weigh it against inertia (no s): the κ needed to
/// reach a target then shrinks with the time step instead of growing like
/// 1/Δt. The multipliers carry s like contact and friction do, so λ is the
/// constraint's reaction force (and λ_A its generalized force on vec(A)) in
/// physical units. Reaction forces change little from one step to the next,
/// so we keep the multipliers across steps: each step starts from the
/// previous step's reaction, and the penalty only has to absorb the change.
class AugmentedLagrangian {
public:
    struct Params {
        double initial_penalty = 1e3;
        double max_penalty = 1e8;
        double satisfied_progress = 0.999;
        double stall_progress = 0.99;
        /// Between the stall and satisfied thresholds, the penalty still
        /// doubles when the remaining gap 1 − η fails to shrink below this
        /// fraction of the previous update's gap.
        double min_gap_reduction = 0.5;
    };

    AugmentedLagrangian() = default;

    explicit AugmentedLagrangian(const Params& params) : m_params(params) { }

    /// @brief Start a new step: capture the start state and targets, reset
    /// the penalties and satisfied flags, and keep each KINEMATIC body's
    /// multipliers from the previous step (others are zeroed).
    /// @param bodies The collection of bodies.
    /// @param x0 Step-start affine DOFs (dim + dim² per body).
    /// @param targets Per-body target rigid poses (only KINEMATIC used).
    /// @param multiplier_scaling The acceleration scaling s = (βΔt)² that
    ///     makes the multipliers physical forces.
    void init(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x0,
        const std::vector<rigid::Pose>& targets,
        const double multiplier_scaling = 1);

    /// @brief Zero one body's multipliers (e.g., when its motion is replaced,
    /// so the previous reaction force no longer applies).
    void reset_multipliers(const size_t body_id);

    /// @brief The constraint's reaction force on a body's position, as of the
    /// last update (zero before then). Only valid after init().
    const VectorMax3d& linear_multiplier(const size_t body_id) const
    {
        return m_lambda[body_id];
    }

    /// @brief Whether any kinematic channel is still unsatisfied.
    bool active() const
    {
        return m_num_kinematic > 0
            && !(m_linear_satisfied && m_angular_satisfied);
    }

    /// @brief Linear progress η_q ∈ (-∞, 1] toward the targets.
    double linear_progress(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x) const;

    /// @brief Angular progress η_A ∈ (-∞, 1] toward the targets.
    double angular_progress(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x) const;

    /// @brief Apply the update policy at the (energy-converged) state x. When
    /// a channel becomes satisfied, we also fold the penalty's force into the
    /// multipliers, so they hold the full reaction force for the next step.
    void update(
        const rigid::RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x);

    bool linear_satisfied() const { return m_linear_satisfied; }
    bool angular_satisfied() const { return m_angular_satisfied; }

    double linear_penalty() const { return m_kappa_q; }
    double angular_penalty() const { return m_kappa_A; }

    // ---- Cumulative functions -----------------------------------------------

    double operator()(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x) const;

    Eigen::VectorXd gradient(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x) const;

    /// @brief The Hessian (block diagonal, constant per step, PSD).
    Eigen::SparseMatrix<double> hessian(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x) const;

private:
    /// @brief Target affine DOFs [p̂; vec(R(θ̂))] per body (dim + dim²).
    VectorMax12d target_dof(const size_t body_id) const;

    /// @brief Method-of-multipliers update: fold the linear penalty's force
    /// at x into λ, i.e., λ ← λ − κ_q m (q − q̂) / s.
    void update_linear_multipliers(
        const rigid::RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x);

    /// @brief λ_A ← λ_A − κ_A W (a − â) / s.
    void update_angular_multipliers(
        const rigid::RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x);

    Params m_params;

    size_t m_num_kinematic = 0;

    /// @brief Scene dimension (set at init).
    int m_dim = 3;

    /// @brief Step-start affine DOFs (η denominators).
    Eigen::VectorXd m_x0;
    /// @brief Target rigid poses.
    std::vector<rigid::Pose> m_targets;
    /// @brief Target rotation matrices (cached).
    std::vector<MatrixMax3d> m_target_rotations;
    /// @brief Diagonal of each kinematic body's angular mass block W (cached).
    std::vector<Eigen::VectorXd> m_angular_weights;

    double m_kappa_q = 1e3;
    double m_kappa_A = 1e3; // NOLINT(readability-identifier-naming)

    /// @brief The acceleration scaling s = (βΔt)² on the multiplier terms.
    double m_multiplier_scaling = 1;

    /// @brief Per-body multipliers (λ linear ∈ R^dim, λ_A angular ∈ R^dim²),
    /// kept across steps as the warm start.
    std::vector<VectorMax3d> m_lambda;
    std::vector<Eigen::VectorXd> m_lambda_A; // NOLINT

    bool m_linear_satisfied = false;
    bool m_angular_satisfied = false;

    /// @brief Each channel's progress at the previous update this step.
    double m_prev_eta_q = -std::numeric_limits<double>::infinity();
    double m_prev_eta_A = -std::numeric_limits<double>::infinity(); // NOLINT
};

} // namespace ipc::affine
