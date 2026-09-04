#pragma once

#include <ipc/dynamics/rigid/rigid_bodies.hpp>

namespace ipc::dynamics {
class ImplicitTimeIntegrator;
}

namespace ipc::affine {

/// @brief Class representing the term ...
class InertialTerm {
public:
    InertialTerm(
        const rigid::RigidBodies& bodies,
        const std::shared_ptr<const dynamics::ImplicitTimeIntegrator>&
            time_integrator);

    /// @brief Update the predicted poses of the rigid bodies.
    /// @param bodies The collection of rigid bodies.
    void update(const rigid::RigidBodies& bodies);

    // ---- Cumulative functions -----------------------------------------------

    /// @brief Compute the total energy for all rigid bodies.
    /// @param bodies The collection of rigid bodies.
    /// @param x The DOFs of the rigid bodies, where the first 3 entries are the positions and the last 3 entries are the rotations.
    /// @return The total energy of the rigid bodies.
    double operator()(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x) const;

    /// @brief Compute the gradient of the total energy for all rigid bodies.
    /// @param bodies The collection of rigid bodies.
    /// @param x The DOFs of the rigid bodies, where the first 3 entries are the positions and the last 3 entries are the rotations.
    /// @return The gradient of the total energy of the rigid bodies.
    Eigen::VectorXd gradient(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x) const;

    /// @brief The (constant, sparse) mass matrix — also the Hessian.
    const Eigen::SparseMatrix<double>& mass_matrix() const { return m_mass; }

private:
    /// @brief (Re)build the block-diagonal ABD mass matrix (zero blocks for
    /// non-DYNAMIC bodies).
    void build_mass_matrix(const rigid::RigidBodies& bodies);

    const std::shared_ptr<const dynamics::ImplicitTimeIntegrator>
        time_integrator;

    /// Mass matrix for the entire system (block diagonal with each block
    /// corresponding to a rigid body)
    Eigen::SparseMatrix<double> m_mass;

    /// Cached predicted poses for the rigid body
    Eigen::VectorXd m_x_hat;
};

} // namespace ipc::affine