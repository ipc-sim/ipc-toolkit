#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/esp/esp_collisions.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <map>

namespace ipc {

// Flag to control parallelism in potential evaluation

class ESPPotential {
public:
    ESPPotential(
        const ESPParameters& _params,
        const bool _use_near_far = true)
        : params(_params)
        , use_near_far(_use_near_far)
    {
    }

    virtual ~ESPPotential() = default;

    // -- Cumulative methods ---------------------------------------------------

    /// @brief Compute the potential for a set of collisions.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param X Degrees of freedom of the collision mesh (e.g., vertices or velocities).
    /// @returns The potential for a set of collisions.
    double operator()(
        const ESPCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> X) const;

    /// @brief Compute the gradient of the potential.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param X Degrees of freedom of the collision mesh (e.g., vertices or velocities).
    /// @returns The gradient of the potential w.r.t. X. This will have a size of |X|.
    Eigen::VectorXd gradient(
        const ESPCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> X) const;

    /// @brief Compute the hessian of the potential.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param X Degrees of freedom of the collision mesh (e.g., vertices or velocities).
    /// @param project_hessian_to_psd Make sure the hessian is positive semi-definite.
    /// @returns The Hessian of the potential w.r.t. X. This will have a size of |X|×|X|.
    virtual Eigen::SparseMatrix<double> hessian(
        const ESPCollisions& collisions,
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
        const ESPCollision& collision,
        Eigen::ConstRef<Eigen::VectorXd> positions) const;

    /// @brief Compute the gradient of the potential for a single collision.
    /// @param collision The collision.
    /// @param positions The collision stencil's positions.
    /// @return The gradient of the potential.
    Eigen::VectorXd gradient(
        const ESPCollision& collision,
        Eigen::ConstRef<Eigen::VectorXd> positions) const;

    /// @brief Compute the hessian of the potential for a single collision.
    /// @param collision The collision.
    /// @param positions The collision stencil's positions.
    /// @return The hessian of the potential.
    Eigen::MatrixXd hessian(
        const ESPCollision& collision,
        Eigen::ConstRef<Eigen::VectorXd> positions,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const;

    using CountMap = std::map<index_t, unsigned>;
    const CountMap& get_edge_evaluation_count() const
    {
        return m_edge_evaluation_count;
    }

    bool get_use_near_far() const { return use_near_far; }

protected:
    /// @brief GCP parameters for collision potential
    ESPParameters params;
    /// @brief Whether to normalize quadrature weights so they sum to 1
    const bool use_near_far;

    mutable CountMap m_edge_evaluation_count;
};

} // namespace ipc
