#pragma once

#include <ipc/potentials/potential.hpp>
#include <ipc/collisions/collisions.hpp>

namespace ipc {

/// @brief Base class for distance-based potentials.
class DistanceBasedPotential : public Potential<Collisions> {
    using Super = Potential<Collisions>;

public:
    DistanceBasedPotential() = default;
    virtual ~DistanceBasedPotential() = default;

    // -- Cumulative methods ---------------------------------------------------

    // NOTE: X in this context are vertex positions.

    using Super::operator();
    using Super::gradient;
    using Super::hessian;

    /// @brief Compute the shape derivative of the potential.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @throws std::runtime_error If the collision collisions were not built with shape derivatives enabled.
    /// @returns The derivative of the force with respect to X, the rest vertices.
    Eigen::SparseMatrix<double> shape_derivative(
        const Collisions& collisions,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices) const;

    // -- Single collision methods ---------------------------------------------

    /// @brief Compute the potential for a single collision.
    /// @param collision The collision.
    /// @param positions The collision stencil's positions.
    /// @return The potential.
    double operator()(const Collision& collision, const VectorMax12d& positions)
        const override;

    /// @brief Compute the gradient of the potential for a single collision.
    /// @param collision The collision.
    /// @param positions The collision stencil's positions.
    /// @return The gradient of the potential.
    VectorMax12d gradient(
        const Collision& collision,
        const VectorMax12d& positions) const override;

    /// @brief Compute the hessian of the potential for a single collision.
    /// @param collision The collision.
    /// @param positions The collision stencil's positions.
    /// @return The hessian of the potential.
    MatrixMax12d hessian(
        const Collision& collision,
        const VectorMax12d& positions,
        const bool project_hessian_to_psd = false) const override;

    /// @brief Compute the shape derivative of the potential for a single collision.
    /// @param[in] collision The collision.
    /// @param[in] vertex_ids The collision stencil's vertex ids.
    /// @param[in] rest_positions The collision stencil's rest positions.
    /// @param[in] positions The collision stencil's positions.
    /// @param[in,out] out Store the triplets of the shape derivative here.
    void shape_derivative(
        const Collision& collision,
        const std::array<long, 4>& vertex_ids,
        const VectorMax12d& rest_positions,
        const VectorMax12d& positions,
        std::vector<Eigen::Triplet<double>>& out) const;

protected:
    /// @brief Compute the unmollified distance-based potential for a collisions.
    /// @param distance_sqr The distance (squared) between the two objects.
    /// @param dmin The minimum distance (unsquared) between the two objects.
    /// @return The unmollified distance-based potential.
    virtual double distance_based_potential(
        const double distance_sqr, const double dmin = 0) const = 0;

    /// @brief Compute the gradient of the unmollified distance-based potential for a collision.
    /// @param distance_sqr The distance (squared) between the two objects.
    /// @param dmin The minimum distance (unsquared) between the two objects.
    /// @return The gradient of the unmollified distance-based potential.
    virtual double distance_based_potential_gradient(
        const double distance_sqr, const double dmin = 0) const = 0;

    /// @brief Compute the hessian of the unmollified distance-based potential for a collision.
    /// @param distance_sqr The distance (squared) between the two objects.
    /// @param dmin The minimum distance (unsquared) between the two objects.
    /// @return The hessian of the unmollified distance-based potential.
    virtual double distance_based_potential_hessian(
        const double distance_sqr, const double dmin = 0) const = 0;
};

} // namespace ipc