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

    using Super::operator();
    using Super::gradient;
    using Super::hessian;

protected:
    /// @brief Compute the potential for a single contact.
    /// @param contact The contact.
    /// @param x The vector of degrees of freedom.
    /// @return The potential.
    double
    operator()(const Contact& contact, const VectorMax12d& x) const override;

    /// @brief Compute the gradient of the potential for a single contact.
    /// @param contact The contact.
    /// @param x The vector of degrees of freedom.
    /// @return The gradient of the potential.
    VectorMax12d
    gradient(const Contact& contact, const VectorMax12d& x) const override;

    /// @brief Compute the hessian of the potential for a single contact.
    /// @param contact The contact.
    /// @param x The vector of degrees of freedom.
    /// @return The hessian of the potential.
    MatrixMax12d hessian(
        const Contact& contact,
        const VectorMax12d& x,
        const bool project_hessian_to_psd = true) const override;

    // ------------------------------------------------------------------------
    // Child classes must implement these functions

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

    // ------------------------------------------------------------------------
};

} // namespace ipc