#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Base class for potentials.
/// @tparam TCollisions The type of the collisions.
template <class TCollisions> class Potential {
protected:
    using TCollision = typename TCollisions::value_type;

public:
    Potential() { }
    virtual ~Potential() { }

    constexpr static int element_size = 3 * TCollision::element_size;

    // -- Cumulative methods ---------------------------------------------------

    /// @brief Compute the potential for a set of collisions.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param X Degrees of freedom of the collision mesh (e.g., vertices or velocities).
    /// @returns The potential for a set of collisions.
    double operator()(
        const TCollisions& collisions,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& X) const;

    /// @brief Compute the gradient of the potential.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param X Degrees of freedom of the collision mesh (e.g., vertices or velocities).
    /// @returns The gradient of the potential w.r.t. X. This will have a size of |X|.
    Eigen::VectorXd gradient(
        const TCollisions& collisions,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& X) const;

    /// @brief Compute the hessian of the potential.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param X Degrees of freedom of the collision mesh (e.g., vertices or velocities).
    /// @param project_hessian_to_psd Make sure the hessian is positive semi-definite.
    /// @returns The Hessian of the potential w.r.t. X. This will have a size of |X|×|X|.
    virtual Eigen::SparseMatrix<double> hessian(
        const TCollisions& collisions,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& X,
        const bool project_hessian_to_psd = false) const;

    // -- Single collision methods ---------------------------------------------

    /// @brief Compute the potential for a single collision.
    /// @param collision The collision.
    /// @param x The collision stencil's degrees of freedom.
    /// @return The potential.
    virtual double
    operator()(const TCollision& collision, const Vector<double, -1, element_size>& x) const = 0;

    /// @brief Compute the gradient of the potential for a single collision.
    /// @param collision The collision.
    /// @param x The collision stencil's degrees of freedom.
    /// @return The gradient of the potential.
    virtual Vector<double, -1, element_size>
    gradient(const TCollision& collision, const Vector<double, -1, element_size>& x) const = 0;

    /// @brief Compute the hessian of the potential for a single collision.
    /// @param collision The collision.
    /// @param x The collision stencil's degrees of freedom.
    /// @return The hessian of the potential.
    virtual MatrixMax<double, element_size, element_size> hessian(
        const TCollision& collision,
        const Vector<double, -1, element_size>& x,
        const bool project_hessian_to_psd = false) const = 0;
};

} // namespace ipc

#include "potential.tpp"