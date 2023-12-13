#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

template <class Contacts> class Potential {
protected:
    using Contact = typename Contacts::value_type;

public:
    Potential() { }
    virtual ~Potential() { }

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
    virtual double
    operator()(const Contact& contact, const VectorMax12d& x) const = 0;

    /// @brief Compute the gradient of the potential for a single contact.
    /// @param contact The contact.
    /// @param x The vector of degrees of freedom.
    /// @return The gradient of the potential.
    virtual VectorMax12d
    gradient(const Contact& contact, const VectorMax12d& x) const = 0;

    /// @brief Compute the hessian of the potential for a single contact.
    /// @param contact The contact.
    /// @param x The vector of degrees of freedom.
    /// @return The hessian of the potential.
    virtual MatrixMax12d hessian(
        const Contact& contact,
        const VectorMax12d& x,
        const bool project_hessian_to_psd = true) const = 0;
};

} // namespace ipc

#include "potential.tpp"