#pragma once

#include <ipc/potentials/distance_based_potential.hpp>

namespace ipc {

/// @brief The normal adhesion potential.
class NormalAdhesionPotential : public DistanceBasedPotential {
public:
    NormalAdhesionPotential(
        const double _dhat_p,
        const double _dhat_a,
        const double _Y,
        const double _eps_c);

    double dhat_p;
    double dhat_a;
    double Y;
    double eps_c;

protected:
    /// @brief Compute the barrier potential for a collision.
    /// @param distance_sqr The distance (squared) between the two objects.
    /// @param dmin The minimum distance (unsquared) between the two objects.
    /// @return The barrier potential.
    double distance_based_potential(
        const double distance_sqr, const double dmin = 0) const override;

    /// @brief Compute the gradient of the barrier potential for a collision.
    /// @param distance_sqr The distance (squared) between the two objects.
    /// @param dmin The minimum distance (unsquared) between the two objects.
    /// @return The gradient of the barrier potential.
    double distance_based_potential_gradient(
        const double distance_sqr, const double dmin = 0) const override;

    /// @brief Compute the hessian of the barrier potential for a collision.
    /// @param distance_sqr The distance (squared) between the two objects.
    /// @param dmin The minimum distance (unsquared) between the two objects.
    /// @return The hessian of the barrier potential.
    double distance_based_potential_hessian(
        const double distance_sqr, const double dmin = 0) const override;
};

} // namespace ipc