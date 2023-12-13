#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/collision_constraints.hpp>

namespace ipc {

using Contacts = CollisionConstraints;
using Contact = CollisionConstraint;

class DistanceBasedPotential {
public:
    DistanceBasedPotential() { }
    virtual ~DistanceBasedPotential() { }

    /// @brief Compute the barrier potential for a set of contacts.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param contacts The set of contacts.
    /// @param dhat The activation distance of the barrier.
    /// @returns The sum of all barrier potentials (not scaled by the barrier stiffness).
    double operator()(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const Contacts& contacts) const;

    /// @brief Compute the gradient of the barrier potential.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param contacts The set of contacts.
    /// @param dhat The activation distance of the barrier.
    /// @returns The gradient of all barrier potentials (not scaled by the barrier stiffness). This will have a size of |vertices|.
    Eigen::VectorXd gradient(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const Contacts& contacts) const;

    /// @brief Compute the hessian of the barrier potential.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param contacts The set of contacts.
    /// @param dhat The activation distance of the barrier.
    /// @param project_hessian_to_psd Make sure the hessian is positive semi-definite.
    /// @returns The hessian of all barrier potentials (not scaled by the barrier stiffness). This will have a size of |vertices|x|vertices|.
    Eigen::SparseMatrix<double> hessian(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const Contacts& contacts,
        const bool project_hessian_to_psd = false) const;

    /// @todo Add shape_derivative() (i.e., the derivative of -gradient() wrt the shape).

protected:
    /// @brief Compute the potential for a single contact.
    /// @param contact The contact.
    /// @param x The vector of degrees of freedom.
    /// @return The potential.
    double potential(const Contact& contact, const VectorMax12d& x) const;

    /// @brief Compute the gradient of the potential for a single contact.
    /// @param contact The contact.
    /// @param x The vector of degrees of freedom.
    /// @return The gradient of the potential.
    VectorMax12d
    potential_gradient(const Contact& contact, const VectorMax12d& x) const;

    /// @brief Compute the hessian of the potential for a single contact.
    /// @param contact The contact.
    /// @param x The vector of degrees of freedom.
    /// @return The hessian of the potential.
    MatrixMax12d potential_hessian(
        const Contact& contact,
        const VectorMax12d& x,
        const bool project_hessian_to_psd = true) const;

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