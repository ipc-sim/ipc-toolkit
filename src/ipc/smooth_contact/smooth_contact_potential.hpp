#pragma once

#include <ipc/potentials/potential.hpp>
#include <ipc/smooth_contact/smooth_collisions.hpp>

namespace ipc {

template <class TCollisions>
class SmoothContactPotential : public Potential<TCollisions> {
    using Super = Potential<TCollisions>;

public:
    SmoothContactPotential(ParameterType& _params) : params(_params) { }

    using Super::element_size;
    using Super::operator();
    using Super::gradient;
    using Super::hessian;

    // -- Single collision methods ---------------------------------------------

    /// @brief Compute the potential for a single collision.
    /// @param collision The collision.
    /// @param positions The collision stencil's positions.
    /// @return The potential.
    double operator()(
        const typename TCollisions::value_type& collision,
        const Vector<double, -1, element_size>& positions) const override
    {
        return collision.weight * collision(positions, params);
    }

    /// @brief Compute the gradient of the potential for a single collision.
    /// @param collision The collision.
    /// @param positions The collision stencil's positions.
    /// @return The gradient of the potential.
    Vector<double, -1, element_size> gradient(
        const typename TCollisions::value_type& collision,
        const Vector<double, -1, element_size>& positions) const override
    {
        return collision.weight * collision.gradient(positions, params);
    }

    /// @brief Compute the hessian of the potential for a single collision.
    /// @param collision The collision.
    /// @param positions The collision stencil's positions.
    /// @return The hessian of the potential.
    MatrixMax<double, element_size, element_size> hessian(
        const typename TCollisions::value_type& collision,
        const Vector<double, -1, element_size>& positions,
        const bool project_hessian_to_psd = false) const override
    {
        MatrixMax<double, element_size, element_size> hess =
            collision.weight * collision.hessian(positions, params);
        return project_hessian_to_psd ? project_to_psd(hess) : hess;
    }

    Eigen::SparseMatrix<double> hessian(
        const TCollisions& collisions,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& X,
        const bool project_hessian_to_psd = false) const override
    {
        return Super::hessian(collisions, mesh, X, project_hessian_to_psd);
    }

protected:
    ParameterType params;
};

} // namespace ipc
