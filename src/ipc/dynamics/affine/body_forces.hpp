#pragma once

#include <ipc/dynamics/rigid/rigid_bodies.hpp>

namespace ipc::dynamics {
class ImplicitTimeIntegrator;
}

namespace ipc::affine {

/// @brief Gravitational and external forces/torques on affine bodies.
///
/// The energy is linear in the affine DOFs: E(x) = wᵀx with a per-body wrench
/// w = [−dt²(m g + f_ext); −dt² vec(τ)] (τ transformed as in the rigid
/// BodyForces, column-major).
class BodyForces {
public:
    explicit BodyForces(
        const std::shared_ptr<const dynamics::ImplicitTimeIntegrator>&
            time_integrator);

    VectorMax3d gravity() const { return m_gravity; }
    void set_gravity(Eigen::ConstRef<VectorMax3d> gravity)
    {
        m_gravity = gravity;
    }

    /// @brief Recompute the per-body wrenches from the current state.
    void update(const rigid::RigidBodies& bodies);

    double operator()(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x) const
    {
        return m_wrench.dot(x);
    }

    Eigen::VectorXd gradient(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x) const
    {
        return m_wrench;
    }

    // The Hessian is identically zero.

private:
    const std::shared_ptr<const dynamics::ImplicitTimeIntegrator>
        time_integrator;

    /// @brief Gravitational acceleration
    VectorMax3d m_gravity = Eigen::Vector3d(0, -9.81, 0);

    /// @brief Stacked per-body wrenches (12 DOF per body in 3D)
    Eigen::VectorXd m_wrench;
};

} // namespace ipc::affine
