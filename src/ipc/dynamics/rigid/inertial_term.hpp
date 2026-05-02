#pragma once

#include <ipc/dynamics/rigid/rigid_potential.hpp>

namespace ipc::rigid {

class ImplicitEuler;

/// @brief Class representing the term ½m‖q - q̂‖² + ½tr((Q - Q̂) J (Q - Q̂)ᵀ)
class InertialTerm : public RigidPotential {
public:
    InertialTerm(const std::shared_ptr<const ImplicitEuler>& _time_integrator)
        : time_integrator(_time_integrator)
    {
    }

    /// @brief Update the predicted poses of the rigid bodies.
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

    // ---- Predicted poses ----------------------------------------------------

    /// @brief Get the predicted poses of the rigid bodies.
    /// @return A vector of predicted poses for each rigid body.
    const std::vector<AffinePose>& predicted_poses() const
    {
        return m_predicted_poses;
    }

private:
    const std::shared_ptr<const ImplicitEuler> time_integrator;

    /// Cached predicted poses for the rigid body
    std::vector<AffinePose> m_predicted_poses;
};

} // namespace ipc::rigid