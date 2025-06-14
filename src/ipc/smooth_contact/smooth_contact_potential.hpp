#pragma once

#include <ipc/potentials/potential.hpp>
#include <ipc/smooth_contact/smooth_collisions.hpp>

namespace ipc {

template <class TCollisions>
class SmoothContactPotential : public Potential<TCollisions> {
    using Super = Potential<TCollisions>;
    using TCollision = typename TCollisions::value_type;

public:
    SmoothContactPotential(const ParameterType& _params) : params(_params) { }

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
        const TCollision& collision,
        Eigen::ConstRef<Vector<double, -1, element_size>> positions)
        const override
    {
        return collision.weight * collision(positions, params);
    }

    /// @brief Compute the gradient of the potential for a single collision.
    /// @param collision The collision.
    /// @param positions The collision stencil's positions.
    /// @return The gradient of the potential.
    Vector<double, -1, element_size> gradient(
        const TCollision& collision,
        Eigen::ConstRef<Vector<double, -1, element_size>> positions)
        const override
    {
        return collision.weight * collision.gradient(positions, params);
    }

    /// @brief Compute the hessian of the potential for a single collision.
    /// @param collision The collision.
    /// @param positions The collision stencil's positions.
    /// @return The hessian of the potential.
    MatrixMax<double, element_size, element_size> hessian(
        const TCollision& collision,
        Eigen::ConstRef<Vector<double, -1, element_size>> positions,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const override
    {
        MatrixMax<double, element_size, element_size> hess =
            collision.weight * collision.hessian(positions, params);
        return project_to_psd(hess, project_hessian_to_psd);
    }

protected:
    ParameterType params;
};

} // namespace ipc
