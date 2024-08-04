#pragma once

#include <ipc/potentials/distance_based_potential.hpp>

namespace ipc {

/// @brief The normal adhesion potential.
class NormalAdhesionPotential : public DistanceBasedPotential {
    using Super = DistanceBasedPotential;

public:
    NormalAdhesionPotential(
        const double _dhat_p,
        const double _dhat_a,
        const double _Y,
        const double _eps_c);

    using Super::operator();
    using Super::gradient;
    using Super::hessian;

    /// @brief The distance of largest adhesion force (\f(\hat{d}_p\f)) (\f(0 < \hat{d}_p < \hat{d}_a\f)).
    double dhat_p;
    /// @brief The adhesion activation distance (\f(\hat{d}_a\f)).
    double dhat_a;
    /// @brief The Young's modulus (\f(Y\f)).
    double Y;
    /// @brief The critical strain (\f(\varepsilon_c\f)).
    double eps_c;

protected:
    /// @brief Compute the barrier potential for a collision.
    /// @param distance_sqr The distance (squared) between the two objects.
    /// @param dmin The minimum distance (unsquared) between the two objects.
    /// @return The barrier potential.
    double
    operator()(const double distance_sqr, const double dmin = 0) const override;

    /// @brief Compute the gradient of the barrier potential for a collision.
    /// @param distance_sqr The distance (squared) between the two objects.
    /// @param dmin The minimum distance (unsquared) between the two objects.
    /// @return The gradient of the barrier potential.
    double
    gradient(const double distance_sqr, const double dmin = 0) const override;

    /// @brief Compute the hessian of the barrier potential for a collision.
    /// @param distance_sqr The distance (squared) between the two objects.
    /// @param dmin The minimum distance (unsquared) between the two objects.
    /// @return The hessian of the barrier potential.
    double
    hessian(const double distance_sqr, const double dmin = 0) const override;
};

} // namespace ipc