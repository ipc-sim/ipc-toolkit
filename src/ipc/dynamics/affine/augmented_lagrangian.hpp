#pragma once

#include <ipc/dynamics/rigid/rigid_bodies.hpp>

namespace ipc::affine {

/// @brief Augmented Lagrangian driving KINEMATIC affine bodies to per-step
/// target poses (the affine analogue of rigid::AugmentedLagrangian).
///
/// With y = [p; vec(A)] per body and the exact ABD mass blocks
/// M = blkdiag(m I₃, W) with W = J ⊗ I₃ (diagonal), the two channels are
/// \f[
///     E = \tfrac{\kappa_q}{2} m \|q - \hat{q}\|^2
///       - \sqrt{m}\, \lambda^\top (q - \hat{q})
///       + \tfrac{\kappa_Q}{2} (a - \hat{a})^\top W (a - \hat{a})
///       - \lambda_A^\top W^{1/2} (a - \hat{a}),
/// \f]
/// where â = vec(R(θ̂)) embeds the rigid target rotation. Fully quadratic
/// with constant PSD Hessian. Same method-of-multipliers schedule as the
/// rigid AL. Pure potential: the simulator applies the (βΔt)² scaling.
class AugmentedLagrangian {
public:
    struct Params {
        double initial_penalty = 1e3;
        double max_penalty = 1e8;
        double satisfied_progress = 0.999;
        double stall_progress = 0.99;
    };

    AugmentedLagrangian() = default;

    explicit AugmentedLagrangian(const Params& params) : m_params(params) { }

    /// @brief Start a new step: capture the start state and targets, reset
    /// the penalties, multipliers, and satisfied flags.
    /// @param bodies The collection of bodies.
    /// @param x0 Step-start affine DOFs (dim + dim² per body).
    /// @param targets Per-body target rigid poses (only KINEMATIC used).
    void init(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x0,
        const std::vector<rigid::Pose>& targets);

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

    /// @brief Apply the update policy at the (energy-converged) state x.
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

    double m_kappa_q = 1e3;
    double m_kappa_A = 1e3; // NOLINT(readability-identifier-naming)

    /// @brief Per-body multipliers (λ linear ∈ R^dim, λ_A angular ∈ R^dim²).
    std::vector<VectorMax3d> m_lambda;
    std::vector<Eigen::VectorXd> m_lambda_A; // NOLINT

    bool m_linear_satisfied = false;
    bool m_angular_satisfied = false;
};

} // namespace ipc::affine
