#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/utils/eigen_ext.hpp>
#include <ipc/smooth_contact/smooth_collisions.hpp>

namespace ipc {

class SmoothContactPotential {
public:
    SmoothContactPotential(const ParameterType& _params) : params(_params) { }
    virtual ~SmoothContactPotential() = default;

    // -- Cumulative methods ---------------------------------------------------

    /// @brief Compute the potential for a set of collisions.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param X Degrees of freedom of the collision mesh (e.g., vertices or velocities).
    /// @returns The potential for a set of collisions.
    double operator()(
        const SmoothCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> X) const;

    /// @brief Compute the gradient of the potential.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param X Degrees of freedom of the collision mesh (e.g., vertices or velocities).
    /// @returns The gradient of the potential w.r.t. X. This will have a size of |X|.
    Eigen::VectorXd gradient(
        const SmoothCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> X) const;

    /// @brief Compute the hessian of the potential.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param X Degrees of freedom of the collision mesh (e.g., vertices or velocities).
    /// @param project_hessian_to_psd Make sure the hessian is positive semi-definite.
    /// @returns The Hessian of the potential w.r.t. X. This will have a size of |X|×|X|.
    virtual Eigen::SparseMatrix<double> hessian(
        const SmoothCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> X,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const;

    // -- Single collision methods ---------------------------------------------

    /// @brief Compute the potential for a single collision.
    /// @param collision The collision.
    /// @param positions The collision stencil's positions.
    /// @return The potential.
    double operator()(
        const SmoothCollision& collision,
        Eigen::ConstRef<Eigen::VectorXd> positions) const;

    /// @brief Compute the gradient of the potential for a single collision.
    /// @param collision The collision.
    /// @param positions The collision stencil's positions.
    /// @return The gradient of the potential.
    Eigen::VectorXd gradient(
        const SmoothCollision& collision,
        Eigen::ConstRef<Eigen::VectorXd> positions) const;

    /// @brief Compute the hessian of the potential for a single collision.
    /// @param collision The collision.
    /// @param positions The collision stencil's positions.
    /// @return The hessian of the potential.
    Eigen::MatrixXd hessian(
        const SmoothCollision& collision,
        Eigen::ConstRef<Eigen::VectorXd> positions,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const;

protected:
    ParameterType params;
};

} // namespace ipc
