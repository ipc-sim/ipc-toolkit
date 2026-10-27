#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/esp/esp_collisions.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class ESPPotential {
public:
    /// @brief Construct an ESP potential.
    /// @param _params ESP parameters.
    /// @param _use_near_far Normalize the quadrature weights (3D only); see
    ///        use_near_far().
    ESPPotential(const ESPParameters& _params, const bool _use_near_far = true)
        : m_params(_params)
        , m_use_near_far(_use_near_far)
    {
        // Automatic lattice rules always weigh their corners.
        const auto& cw = _params.corner_weights();
        const bool corners = _params.lattice_spacing() > 0
            || (cw[0] > 0 && cw[1] > 0 && cw[2] > 0);
        if (_params.local_quadrature_transfer && (!_use_near_far || !corners)) {
            throw std::invalid_argument(
                "Local quadrature transfer requires normalized near/far "
                "integration and a face rule with positive weight at every "
                "corner.");
        }
    }

    virtual ~ESPPotential() = default;

    // -- Cumulative methods ---------------------------------------------------

    /// @brief Compute the potential for a set of collisions.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param X Degrees of freedom of the collision mesh (e.g., vertices or velocities).
    /// @returns The potential for a set of collisions.
    /// @throws std::invalid_argument if the potential's quadrature differs
    ///         from the one collisions were built with (see
    ///         ESPCollisions::check_quadrature()).
    double operator()(
        const ESPCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> X) const;

    /// @brief Compute the gradient of the potential.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param X Degrees of freedom of the collision mesh (e.g., vertices or velocities).
    /// @returns The gradient of the potential w.r.t. X. This will have a size of |X|.
    /// @throws std::invalid_argument if the potential's quadrature differs
    ///         from the one collisions were built with (see
    ///         ESPCollisions::check_quadrature()).
    Eigen::VectorXd gradient(
        const ESPCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> X) const;

    /// @brief Compute the hessian of the potential.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param X Degrees of freedom of the collision mesh (e.g., vertices or velocities).
    /// @param project_hessian_to_psd Make sure the hessian is positive semi-definite.
    ///        With normalized weights and edge-edge terms, each face's block is
    ///        projected as a whole, since its weight-derivative terms are
    ///        indefinite; otherwise each stencil's block is projected.
    /// @returns The Hessian of the potential w.r.t. X. This will have a size of |X|×|X|.
    /// @throws std::invalid_argument if the potential's quadrature differs
    ///         from the one collisions were built with (see
    ///         ESPCollisions::check_quadrature()).
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
    /// @param project_hessian_to_psd Make sure the hessian is positive semi-definite.
    /// @return The hessian of the potential.
    Eigen::MatrixXd hessian(
        const ESPCollision& collision,
        Eigen::ConstRef<Eigen::VectorXd> positions,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const;

    /// @brief Whether to normalize the quadrature weights so they sum to 1,
    ///        splitting near and far contributions when 0 < dbar_factor <= 2.
    ///        The far part is nonzero for dbar_factor < 2 (e.g., at 1) and
    ///        vanishes identically at 2, where the split equals the unsplit
    ///        normalized potential. 3D only: the 2D potential is never
    ///        normalized.
    bool use_near_far() const { return m_use_near_far; }

protected:
    /// @brief ESP parameters for collision potential
    ESPParameters m_params;
    /// @brief See use_near_far().
    const bool m_use_near_far;
};

} // namespace ipc
