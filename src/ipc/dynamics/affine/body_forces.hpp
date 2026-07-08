#pragma once

#include <ipc/dynamics/rigid/rigid_bodies.hpp>

namespace ipc::affine {

/// @brief Gravitational and external forces/torques on affine bodies.
///
/// The energy is linear in the affine DOFs: E(x) = wᵀx with a per-body wrench
/// w = [−(m g + f_ext); −vec(τ)] (τ transformed as in the rigid BodyForces,
/// column-major).
///
/// This is a pure potential (no time-integrator scaling): the simulator is
/// responsible for scaling it into the incremental potential (by (βΔt)²).
class BodyForces {
public:
    BodyForces() = default;

    VectorMax3d gravity() const { return m_gravity; }
    void set_gravity(Eigen::ConstRef<VectorMax3d> gravity)
    {
        m_gravity = gravity;
    }

    /// @brief Recompute the per-body wrenches from the current state.
    /// @param bodies The collection of bodies.
    /// @param poses The poses of the bodies at the start of the step (used to
    ///     transform world-space torques into body space).
    void update(
        const rigid::RigidBodies& bodies,
        const std::vector<rigid::AffinePose>& poses);

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
    /// @brief Gravitational acceleration
    VectorMax3d m_gravity = Eigen::Vector3d(0, -9.81, 0);

    /// @brief Stacked per-body wrenches (dim + dim² DOF per body)
    Eigen::VectorXd m_wrench;
};

} // namespace ipc::affine
