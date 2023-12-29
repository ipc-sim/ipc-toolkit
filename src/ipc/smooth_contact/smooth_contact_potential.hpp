#pragma once

#include <ipc/potentials/potential.hpp>
#include <ipc/collisions/collisions.hpp>

namespace ipc {

class SmoothContactPotential : public Potential<Collisions> {
    using Super = Potential<Collisions>;

public:
    SmoothContactPotential(ParameterType &_params) : params(_params) { }
    virtual ~SmoothContactPotential() { }

    // -- Cumulative methods ---------------------------------------------------

    // NOTE: X in this context are vertex positions.

    using Super::operator();
    using Super::gradient;
    using Super::hessian;

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

protected:
    ParameterType params;

};

} // namespace ipc