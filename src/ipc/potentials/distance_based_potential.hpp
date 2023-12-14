#pragma once

#include <ipc/potentials/potential.hpp>
#include <ipc/collisions/collision_constraints.hpp>

namespace ipc {

class DistanceBasedPotential : public Potential<CollisionConstraints> {
    using Super = Potential<CollisionConstraints>;
    using Contacts = CollisionConstraints;
    using Contact = Super::Contact;

public:
    DistanceBasedPotential() { }
    virtual ~DistanceBasedPotential() { }

    // -- Cumulative methods ---------------------------------------------------

    // NOTE: X in this context are vertex positions.

    using Super::operator();
    using Super::gradient;
    using Super::hessian;

    /// @brief Compute the shape derivative of the potential.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param contacts The set of contacts.
    /// @throws std::runtime_error If the collision constraints were not built with shape derivatives enabled.
    /// @returns The derivative of the force with respect to X, the rest vertices.
    Eigen::SparseMatrix<double> shape_derivative(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const Contacts& contacts) const;

    // -- Single contact methods -----------------------------------------------

    /// @brief Compute the potential for a single contact.
    /// @param contact The contact.
    /// @param positions The contact stencil's positions.
    /// @return The potential.
    double operator()(
        const Contact& contact, const VectorMax12d& positions) const override;

    /// @brief Compute the gradient of the potential for a single contact.
    /// @param contact The contact.
    /// @param positions The contact stencil's positions.
    /// @return The gradient of the potential.
    VectorMax12d gradient(
        const Contact& contact, const VectorMax12d& positions) const override;

    /// @brief Compute the hessian of the potential for a single contact.
    /// @param contact The contact.
    /// @param positions The contact stencil's positions.
    /// @return The hessian of the potential.
    MatrixMax12d hessian(
        const Contact& contact,
        const VectorMax12d& positions,
        const bool project_hessian_to_psd = true) const override;

    /// @brief Compute the shape derivative of the potential for a single contact.
    /// @param[in] contact The contact.
    /// @param[in] vertex_ids The contact stencil's vertex ids.
    /// @param[in] rest_positions The contact stencil's rest positions.
    /// @param[in] positions The contact stencil's positions.
    /// @param[in,out] out Store the triplets of the shape derivative here.
    void shape_derivative(
        const Contact& contact,
        const std::array<long, 4>& vertex_ids,
        const VectorMax12d& rest_positions,
        const VectorMax12d& positions,
        std::vector<Eigen::Triplet<double>>& out) const;

protected:
    /// @brief Compute the unmollified distance-based potential for a contacts.
    /// @param distance_sqr The distance (squared) between the two objects.
    /// @param dmin The minimum distance (unsquared) between the two objects.
    /// @return The unmollified distance-based potential.
    virtual double distance_based_potential(
        const double distance_sqr, const double dmin = 0) const = 0;

    /// @brief Compute the gradient of the unmollified distance-based potential for a contact.
    /// @param distance_sqr The distance (squared) between the two objects.
    /// @param dmin The minimum distance (unsquared) between the two objects.
    /// @return The gradient of the unmollified distance-based potential.
    virtual double distance_based_potential_gradient(
        const double distance_sqr, const double dmin = 0) const = 0;

    /// @brief Compute the hessian of the unmollified distance-based potential for a contact.
    /// @param distance_sqr The distance (squared) between the two objects.
    /// @param dmin The minimum distance (unsquared) between the two objects.
    /// @return The hessian of the unmollified distance-based potential.
    virtual double distance_based_potential_hessian(
        const double distance_sqr, const double dmin = 0) const = 0;
};

} // namespace ipc