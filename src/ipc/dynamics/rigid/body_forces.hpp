#pragma once

#include <ipc/dynamics/rigid/pose.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/dynamics/rigid/time_integrator.hpp>

namespace ipc::rigid {

/// @brief Class representing the term qᵀf + tr(Qᵀτ)
class BodyForces {
public:
    BodyForces(const std::shared_ptr<const ImplicitEuler>& _time_integrator)
        : time_integrator(_time_integrator)
    {
    }

    /// @brief Update the forces and torques of the rigid bodies.
    /// @param bodies The collection of rigid bodies.
    void update(const RigidBodies& bodies);

    // ---- Cumulative functions -----------------------------------------------

    /// @brief Compute the total energy for all rigid bodies.
    /// @param bodies The collection of rigid bodies.
    /// @param x The DOFs of the rigid bodies, where the first 3 entries are the positions and the last 3 entries are the rotations.
    /// @return The total energy of the rigid bodies.
    double operator()(
        const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x) const;

    /// @brief Compute the gradient of the total energy for all rigid bodies.
    /// @param bodies The collection of rigid bodies.
    /// @param x The DOFs of the rigid bodies, where the first 3 entries are the positions and the last 3 entries are the rotations.
    /// @return The gradient of the total energy of the rigid bodies.
    Eigen::VectorXd gradient(
        const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x) const;

    /// @brief Compute the Hessian of the total energy for all rigid bodies.
    /// @param bodies The collection of rigid bodies.
    /// @param x The DOFs of the rigid bodies, where the first 3 entries are the positions and the last 3 entries are the rotations.
    /// @return The Hessian of the total energy of the rigid bodies.
    Eigen::MatrixXd hessian(
        const RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const;

    // ---- Per-body functions -------------------------------------------------

    /// @brief Compute the energy of a rigid body at a given pose.
    /// @param body The rigid body.
    /// @param pose The pose of the rigid body.
    /// @param force The linear force on the rigid body.
    /// @param torque The torque on the rigid body.
    /// @return The energy of the rigid body at the given pose.
    double operator()(
        const RigidBody& body,
        Eigen::ConstRef<VectorMax6d> x,
        Eigen::ConstRef<VectorMax3d> force,
        Eigen::ConstRef<MatrixMax3d> torque) const;

    /// @brief Compute the gradient of the energy of a rigid body at a given pose.
    /// @param body The rigid body.
    /// @param x The pose of the rigid body.
    /// @param force The linear force on the rigid body.
    /// @param torque The torque on the rigid body.
    /// @return The gradient of the energy of the rigid body at the given pose.
    VectorMax6d gradient(
        const RigidBody& body,
        Eigen::ConstRef<VectorMax6d> x,
        Eigen::ConstRef<VectorMax3d> force,
        Eigen::ConstRef<MatrixMax3d> torque) const;

    /// @brief Compute the Hessian of the energy of a rigid body at a given pose.
    /// @param body The rigid body.
    /// @param x The pose of the rigid body.
    /// @param force The linear force on the rigid body.
    /// @param torque The torque on the rigid body.
    /// @return The Hessian of the energy of the rigid body at the given pose.
    MatrixMax6d hessian(
        const RigidBody& body,
        Eigen::ConstRef<VectorMax6d> x,
        Eigen::ConstRef<VectorMax3d> force,
        Eigen::ConstRef<MatrixMax3d> torque,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const;

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