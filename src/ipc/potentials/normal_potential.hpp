#pragma once

#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/potentials/potential.hpp>

namespace ipc {

/// @brief Base class for distance-based potentials.
class NormalPotential : public Potential<NormalCollisions> {
    using Super = Potential<NormalCollisions>;

public:
    NormalPotential() = default;
    virtual ~NormalPotential() = default;

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
        const NormalCollisions& collisions,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices) const;

    // -- Single collision methods ---------------------------------------------

    /// @brief Compute the potential for a single collision.
    /// @param collision The collision.
    /// @param positions The collision stencil's positions.
    /// @return The potential.
    double operator()(
        const NormalCollision& collision,
        const VectorMax12d& positions) const override;

    /// @brief Compute the gradient of the potential for a single collision.
    /// @param collision The collision.
    /// @param positions The collision stencil's positions.
    /// @return The gradient of the potential.
    VectorMax12d gradient(
        const NormalCollision& collision,
        const VectorMax12d& positions) const override;

    /// @brief Compute the hessian of the potential for a single collision.
    /// @param collision The collision.
    /// @param positions The collision stencil's positions.
    /// @return The hessian of the potential.
    MatrixMax12d hessian(
        const NormalCollision& collision,
        const VectorMax12d& positions,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const override;

    /// @brief Compute the shape derivative of the potential for a single collision.
    /// @param[in] collision The collision.
    /// @param[in] vertex_ids The collision stencil's vertex ids.
    /// @param[in] rest_positions The collision stencil's rest positions.
    /// @param[in] positions The collision stencil's positions.
    /// @param[in,out] out Store the triplets of the shape derivative here.
    void shape_derivative(
        const NormalCollision& collision,
        const std::array<long, 4>& vertex_ids,
        const VectorMax12d& rest_positions,
        const VectorMax12d& positions,
        std::vector<Eigen::Triplet<double>>& out) const;

    /// @brief Compute the force magnitude for a collision.
    /// @param distance_squared The squared distance between elements.
    /// @param dmin The minimum distance offset to the barrier.
    /// @param barrier_stiffness The barrier stiffness.
    /// @return The force magnitude.
    virtual double force_magnitude(
        const double distance_squared,
        const double dmin,
        const double barrier_stiffness) const = 0;

    /// @brief Compute the gradient of the force magnitude for a collision.
    /// @param distance_squared The squared distance between elements.
    /// @param distance_squared_gradient The gradient of the squared distance.
    /// @param dmin The minimum distance offset to the barrier.
    /// @param barrier_stiffness The stiffness of the barrier.
    /// @return The gradient of the force.
    virtual VectorMax12d force_magnitude_gradient(
        const double distance_squared,
        const VectorMax12d& distance_squared_gradient,
        const double dmin,
        const double barrier_stiffness) const = 0;

protected:
    /// @brief Compute the unmollified distance-based potential for a collisions.
    /// @param distance_squared The distance (squared) between the two objects.
    /// @param dmin The minimum distance (unsquared) between the two objects.
    /// @return The unmollified distance-based potential.
    virtual double
    operator()(const double distance_squared, const double dmin = 0) const = 0;

    /// @brief Compute the gradient of the unmollified distance-based potential for a collision.
    /// @param distance_squared The distance (squared) between the two objects.
    /// @param dmin The minimum distance (unsquared) between the two objects.
    /// @return The gradient of the unmollified distance-based potential.
    virtual double
    gradient(const double distance_squared, const double dmin = 0) const = 0;

    /// @brief Compute the hessian of the unmollified distance-based potential for a collision.
    /// @param distance_squared The distance (squared) between the two objects.
    /// @param dmin The minimum distance (unsquared) between the two objects.
    /// @return The hessian of the unmollified distance-based potential.
    virtual double
    hessian(const double distance_squared, const double dmin = 0) const = 0;
};

} // namespace ipc