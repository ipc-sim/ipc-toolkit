#pragma once

#include <ipc/dynamics/rigid/pose.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>

namespace ipc::rigid {

/// @brief Class representing the potential energy of a rigid body, which includes the inertial term and the external forces.
class RigidPotential {
public:
    virtual ~RigidPotential() = default;

    /// @brief Update the predicted poses of the rigid bodies.
    /// @param bodies The collection of rigid bodies.
    virtual void update(const RigidBodies& bodies) = 0;

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
    /// @param body_id The index of the rigid body in the collection of rigid bodies.
    /// @param body The rigid body.
    /// @param pose The pose of the rigid body.
    /// @return The energy of the rigid body at the given pose.
    virtual double operator()(
        const size_t body_id,
        const RigidBody& body,
        Eigen::ConstRef<VectorMax6d> x) const = 0;

    /// @brief Compute the gradient of the energy of a rigid body at a given pose.
    /// @param body_id The index of the rigid body in the collection of rigid bodies.
    /// @param body The rigid body.
    /// @param x The pose of the rigid body.
    /// @return The gradient of the energy of the rigid body at the given pose.
    virtual VectorMax6d gradient(
        const size_t body_id,
        const RigidBody& body,
        Eigen::ConstRef<VectorMax6d> x) const = 0;

    /// @brief Compute the Hessian of the energy of a rigid body at a given pose.
    /// @param body_id The index of the rigid body in the collection of rigid bodies.
    /// @param body The rigid body.
    /// @param x The pose of the rigid body.
    /// @return The Hessian of the energy of the rigid body at the given pose.
    virtual MatrixMax6d hessian(
        const size_t body_id,
        const RigidBody& body,
        Eigen::ConstRef<VectorMax6d> x,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const = 0;
};

} // namespace ipc::rigid