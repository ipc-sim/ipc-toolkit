#pragma once

#include <ipc/dynamics/rigid/rigid_potential.hpp>

namespace ipc::rigid {

class ImplicitEuler;

/// @brief Class representing the term qᵀf + tr(Qᵀτ)
class BodyForces : public RigidPotential {
public:
    BodyForces(const std::shared_ptr<const ImplicitEuler>& _time_integrator)
        : time_integrator(_time_integrator)
    {
    }

    /// @brief Update the forces and torques of the rigid bodies.
    /// @param bodies The collection of rigid bodies.
    void update(const RigidBodies& bodies) override;

    // ---- Cumulative functions -----------------------------------------------

    using RigidPotential::operator();
    using RigidPotential::gradient;
    using RigidPotential::hessian;

    // ---- Per-body functions -------------------------------------------------

    /// @brief Compute the energy of a rigid body at a given pose.
    /// @param body_id The index of the rigid body in the collection of rigid bodies.
    /// @param body The rigid body.
    /// @param pose The pose of the rigid body.
    /// @return The energy of the rigid body at the given pose.
    double operator()(
        const size_t body_id,
        const RigidBody& body,
        Eigen::ConstRef<VectorMax6d> x) const override;

    /// @brief Compute the gradient of the energy of a rigid body at a given pose.
    /// @param body_id The index of the rigid body in the collection of rigid bodies.
    /// @param body The rigid body.
    /// @param x The pose of the rigid body.
    /// @return The gradient of the energy of the rigid body at the given pose.
    VectorMax6d gradient(
        const size_t body_id,
        const RigidBody& body,
        Eigen::ConstRef<VectorMax6d> x) const override;

    /// @brief Compute the Hessian of the energy of a rigid body at a given pose.
    /// @param body_id The index of the rigid body in the collection of rigid bodies.
    /// @param body The rigid body.
    /// @param x The pose of the rigid body.
    /// @return The Hessian of the energy of the rigid body at the given pose.
    MatrixMax6d hessian(
        const size_t body_id,
        const RigidBody& body,
        Eigen::ConstRef<VectorMax6d> x,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const override;

    // ---- Gravity ------------------------------------------------------------

    const VectorMax3d& gravity() const { return m_gravity; }

    void set_gravity(Eigen::ConstRef<VectorMax3d> gravity)
    {
        m_gravity = gravity;
    }

    const std::vector<VectorMax3d>& forces() const { return m_forces; }
    const std::vector<MatrixMax3d>& torques() const { return m_torques; }

private:
    const std::shared_ptr<const ImplicitEuler> time_integrator;

    std::vector<VectorMax3d> m_forces;
    std::vector<MatrixMax3d> m_torques;

    VectorMax3d m_gravity = VectorMax3d::Zero(3);
};

} // namespace ipc::rigid