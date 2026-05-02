#pragma once

#include <ipc/dynamics/rigid/rigid_bodies.hpp>

namespace ipc::rigid {
class ImplicitEuler;
}

namespace ipc::affine {

/// @brief Class representing the term ...
class InertialTerm {
public:
    InertialTerm(
        const rigid::RigidBodies& bodies,
        const std::shared_ptr<const rigid::ImplicitEuler>& time_integrator);

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

    /// @brief Compute the Hessian of the total energy for all rigid bodies.
    /// @param bodies The collection of rigid bodies.
    /// @param x The DOFs of the rigid bodies, where the first 3 entries are the positions and the last 3 entries are the rotations.
    /// @return The Hessian of the total energy of the rigid bodies.
    Eigen::MatrixXd hessian(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const;

private:
    static Eigen::SparseMatrix<double>
    J(Eigen::ConstRef<Eigen::MatrixXd> rest_positions) // NOLINT
    {
        std::vector<Eigen::Triplet<double>> triplets;

        for (int i = 0; i < rest_positions.rows(); i++) {
            for (int j = 0; j < rest_positions.cols(); j++) {
                triplets.emplace_back(rest_positions.rows() * j + i, j, 1);
            }
        }

        // I ⊗ x̄
        for (int i = 0; i < rest_positions.rows(); i++) {
            for (int j = 0; j < rest_positions.cols(); j++) {
                for (int k = 0; k < 3; k++) {
                    triplets.emplace_back(
                        i + k * rest_positions.rows(),
                        j + k * rest_positions.cols() + 3,
                        rest_positions(i, j));
                }
            }
        }

        Eigen::SparseMatrix<double> J(rest_positions.size(), 12);
        J.setFromTriplets(triplets.begin(), triplets.end());
        return J;
    }

    const std::shared_ptr<const rigid::ImplicitEuler> time_integrator;

    /// Mass matrix for the entire system (block diagonal with each block
    /// corresponding to a rigid body)
    Eigen::SparseMatrix<double> m_mass;

    /// Cached predicted poses for the rigid body
    Eigen::VectorXd m_x_hat;
};

} // namespace ipc::affine