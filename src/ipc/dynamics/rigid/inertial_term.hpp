#pragma once

#include <ipc/dynamics/rigid/pose.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/dynamics/rigid/time_integrator.hpp>

namespace ipc::rigid {

/// @brief Class representing the term ½m‖q - q̂‖² + ½tr((Q - Q̂) J (Q - Q̂)ᵀ)
class InertialTerm {
public:
    InertialTerm(const std::shared_ptr<const ImplicitEuler>& _time_integrator)
        : time_integrator(_time_integrator)
    {
    }

    /// @brief Update the predicted poses of the rigid bodies.
    /// @param bodies The collection of rigid bodies.
    void update(const RigidBodies& bodies);

    // ---- Cumulative functions -----------------------------------------------

    /// @brief Compute the total energy for all rigid bodies.
    /// @param bodies The collection of rigid bodies.
    /// @param x The DOFs of the rigid bodies, where the first 3 entries are the positions and the last 3 entries are the rotations.
    /// @return The total energy of the rigid bodies.
    double
    operator()(const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x);

    /// @brief Compute the gradient of the total energy for all rigid bodies.
    /// @param bodies The collection of rigid bodies.
    /// @param x The DOFs of the rigid bodies, where the first 3 entries are the positions and the last 3 entries are the rotations.
    /// @return The gradient of the total energy of the rigid bodies.
    Eigen::VectorXd
    gradient(const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x);

    /// @brief Compute the Hessian of the total energy for all rigid bodies.
    /// @param bodies The collection of rigid bodies.
    /// @param x The DOFs of the rigid bodies, where the first 3 entries are the positions and the last 3 entries are the rotations.
    /// @return The Hessian of the total energy of the rigid bodies.
    Eigen::MatrixXd
    hessian(const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x);

    // ---- Per-body functions -------------------------------------------------

    /// @brief Compute the energy of a rigid body at a given pose.
    /// @param body The rigid body.
    /// @param pose The pose of the rigid body.
    /// @return The energy of the rigid body at the given pose.
    double operator()(
        const RigidBody& body,
        Eigen::ConstRef<VectorMax6d> x,
        Eigen::ConstRef<VectorMax3d> q_hat,
        Eigen::ConstRef<MatrixMax3d> Q_hat) const;

    /// @brief Compute the gradient of the energy of a rigid body at a given pose.
    /// @param body The rigid body.
    /// @param x The pose of the rigid body.
    /// @param q_hat The predicted rotation of the rigid body.
    /// @param Q_hat The predicted rotation matrix of the rigid body.
    /// @return The gradient of the energy of the rigid body at the given pose.
    VectorMax6d gradient(
        const RigidBody& body,
        Eigen::ConstRef<VectorMax6d> x,
        Eigen::ConstRef<VectorMax3d> q_hat,
        Eigen::ConstRef<MatrixMax3d> Q_hat) const;

    /// @brief Compute the Hessian of the energy of a rigid body at a given pose.
    /// @param body The rigid body.
    /// @param x The pose of the rigid body.
    /// @param q_hat The predicted rotation of the rigid body.
    /// @param Q_hat The predicted rotation matrix of the rigid body.
    /// @return The Hessian of the energy of the rigid body at the given pose.
    MatrixMax6d hessian(
        const RigidBody& body,
        Eigen::ConstRef<VectorMax6d> x,
        Eigen::ConstRef<VectorMax3d> q_hat,
        Eigen::ConstRef<MatrixMax3d> Q_hat) const;

private:
    const std::shared_ptr<const ImplicitEuler> time_integrator;

    /// Cached predicted poses for the rigid body
    std::vector<AffinePose> predicted_poses;
};

} // namespace ipc::rigid