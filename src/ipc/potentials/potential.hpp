#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Base class for potentials.
/// @tparam TCollisions The type of the collisions.
template <class TCollisions> class Potential {
protected:
    using TCollision = typename TCollisions::value_type;
    /// @brief Maximum degrees of freedom per collision
    static constexpr int ELEMENT_SIZE = 3 * TCollision::ELEMENT_SIZE;
    using VectorMaxNd = Vector<double, Eigen::Dynamic, ELEMENT_SIZE>;
    using MatrixMaxNd = MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>;

public:
    Potential() = default;
    virtual ~Potential() = default;

    // -- Cumulative methods ---------------------------------------------------

    /// @brief Compute the potential for a set of collisions.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param X Degrees of freedom of the collision mesh (e.g., vertices or velocities).
    /// @returns The potential for a set of collisions.
    double operator()(
        const TCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> X) const;

    /// @brief Compute the gradient of the potential.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param X Degrees of freedom of the collision mesh (e.g., vertices or velocities).
    /// @returns The gradient of the potential w.r.t. X. This will have a size of |X|.
    Eigen::VectorXd gradient(
        const TCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> X) const;

    /// @brief Compute the hessian of the potential.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param X Degrees of freedom of the collision mesh (e.g., vertices or velocities).
    /// @param project_hessian_to_psd Make sure the hessian is positive semi-definite.
    /// @returns The Hessian of the potential w.r.t. X. This will have a size of |X|×|X|.
    virtual Eigen::SparseMatrix<double> hessian(
        const TCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> X,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const;

    // -- Single collision methods ---------------------------------------------

    /// @brief Compute the potential for a single collision.
    /// @param collision The collision.
    /// @param x The collision stencil's degrees of freedom.
    /// @return The potential.
    virtual double operator()(
        const TCollision& collision, Eigen::ConstRef<VectorMaxNd> x) const = 0;

    /// @brief Compute the gradient of the potential for a single collision.
    /// @param collision The collision.
    /// @param x The collision stencil's degrees of freedom.
    /// @return The gradient of the potential.
    virtual VectorMaxNd gradient(
        const TCollision& collision, Eigen::ConstRef<VectorMaxNd> x) const = 0;

    /// @brief Compute the hessian of the potential for a single collision.
    /// @param collision The collision.
    /// @param x The collision stencil's degrees of freedom.
    /// @param project_hessian_to_psd Whether to project the hessian to the positive semi-definite cone.
    /// @return The hessian of the potential.
    virtual MatrixMaxNd hessian(
        const TCollision& collision,
        Eigen::ConstRef<VectorMaxNd> x,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const = 0;
};

} // namespace ipc