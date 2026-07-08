#pragma once

#include <ipc/dynamics/affine/augmented_lagrangian.hpp>
#include <ipc/dynamics/affine/body_forces.hpp>
#include <ipc/dynamics/affine/inertial_term.hpp>
#include <ipc/dynamics/affine/orthogonality_potential.hpp>
#include <ipc/dynamics/to_affine.hpp>

#include <memory>
#include <optional>

namespace ipc::dynamics {

class ImplicitTimeIntegrator;

/// @brief The body (non-contact) part of the incremental potential, expressed
/// in the optimization DOFs x through a single to-affine map.
///
/// It owns the per-body affine-coordinate terms — inertia, body forces, the
/// SO(dim) orthogonality penalty (affine bodies only), and the augmented
/// Lagrangian — evaluates each at y = to_affine(x), sums their
/// gradients/Hessians (with the incremental-potential scaling: inertia
/// unscaled, the rest × s), and applies the to-affine chain rule once:
///     ∇ₓE  = (dy/dx)ᵀ (g_inertia + s (g_bf + g_orth [+ g_al]))
///     ∇²ₓE = (dy/dx)ᵀ H (dy/dx) + Σₖ gₖ d²yₖ/dx².
///
/// These are the block-diagonal-per-body terms; pairwise contact/friction use
/// the vertex-space chain rule and are added separately.
///
/// The augmented Lagrangian is separable via @p include_al so the adaptive
/// barrier stiffness can be seeded from the physical forces only.
class BodyPotentials {
public:
    BodyPotentials(
        const rigid::RigidBodies& bodies,
        std::shared_ptr<const ImplicitTimeIntegrator> time_integrator,
        std::shared_ptr<const ToAffine> to_affine,
        const double orthogonality_stiffness,
        const affine::AugmentedLagrangian::Params& al_params = {});

    VectorMax3d gravity() const { return m_body_forces.gravity(); }
    void set_gravity(Eigen::ConstRef<VectorMax3d> gravity)
    {
        m_body_forces.set_gravity(gravity);
    }

    /// @brief Refresh the predicted state and body-force torques for a new step.
    void update(const rigid::RigidBodies& bodies);

    // ---- Incremental-potential contribution ---------------------------------

    double energy(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x,
        const double scaling,
        const bool include_al = true) const;

    Eigen::VectorXd gradient(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x,
        const double scaling,
        const bool include_al = true) const;

    Eigen::SparseMatrix<double> hessian(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x,
        const double scaling,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE,
        const bool include_al = true) const;

    // ---- Augmented Lagrangian lifecycle (operates in affine coordinates) ----

    void init_augmented_lagrangian(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x0,
        const std::vector<rigid::Pose>& targets);

    void update_augmented_lagrangian(
        const rigid::RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x);

    bool augmented_lagrangian_active() const { return m_al.active(); }
    bool augmented_lagrangian_linear_satisfied() const
    {
        return m_al.linear_satisfied();
    }
    bool augmented_lagrangian_angular_satisfied() const
    {
        return m_al.angular_satisfied();
    }

    double augmented_lagrangian_linear_progress(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x) const
    {
        return m_al.linear_progress(bodies, m_to_affine->to_affine(x));
    }
    double augmented_lagrangian_angular_progress(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x) const
    {
        return m_al.angular_progress(bodies, m_to_affine->to_affine(x));
    }

    affine::AugmentedLagrangian& augmented_lagrangian() { return m_al; }
    const affine::AugmentedLagrangian& augmented_lagrangian() const
    {
        return m_al;
    }

    // ---- Accessors ----------------------------------------------------------

    const ToAffine& to_affine() const { return *m_to_affine; }

    /// @brief The constant affine-coordinate inertial (mass) matrix.
    const Eigen::SparseMatrix<double>& mass_matrix() const
    {
        return m_inertial.mass_matrix();
    }

private:
    /// @brief The summed affine-coordinate gradient (inertia + s·(rest)).
    Eigen::VectorXd affine_gradient(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> y,
        const double scaling,
        const bool include_al) const;

    std::shared_ptr<const ImplicitTimeIntegrator> m_time_integrator;
    std::shared_ptr<const ToAffine> m_to_affine;

    affine::InertialTerm m_inertial;
    affine::BodyForces m_body_forces;
    affine::AugmentedLagrangian m_al;
    std::optional<affine::OrthogonalityPotential> m_orthogonality;
};

} // namespace ipc::dynamics
