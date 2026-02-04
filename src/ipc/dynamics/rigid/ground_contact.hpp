#pragma once

#include <ipc/barrier/barrier.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/SparseCore>

namespace ipc::rigid {

struct Pose;
class RigidBodies;

class GroundContact {
public:
    GroundContact(const double height);

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
    /// @param project_hessian_to_psd The method to project the Hessian to be positive semi-definite.
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
    /// @return The energy of the rigid body at the given pose.
    double operator()(
        const RigidBodies& bodies,
        const size_t body_index,
        const Pose& pose) const;

    /// @brief Compute the gradient of the energy of a rigid body at a given pose.
    /// @param body The rigid body.
    /// @param x The pose of the rigid body.
    /// @return The gradient of the energy of the rigid body at the given pose.
    VectorMax6d gradient(
        const RigidBodies& bodies,
        const size_t body_index,
        const Pose& pose) const;

    /// @brief Compute the Hessian of the energy of a rigid body at a given pose.
    /// @param body The rigid body.
    /// @param x The pose of the rigid body.
    /// @param project_hessian_to_psd The method to project the Hessian to be positive semi-definite.
    /// @return The Hessian of the energy of the rigid body at the given pose.
    MatrixMax6d hessian(
        const RigidBodies& bodies,
        const size_t body_index,
        const Pose& pose,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const;

    // ---- Gravity ------------------------------------------------------------

    double ground_height() const { return m_ground_height; }

    void set_ground_height(const double ground_height)
    {
        m_ground_height = ground_height;
    }

    double dhat() const { return m_dhat; }

    void set_dhat(const double dhat) { m_dhat = dhat; }

private:
    double m_ground_height;
    std::shared_ptr<Barrier> m_barrier;
    double m_dhat; // Barrier function parameter
};

} // namespace ipc::rigid