#pragma once

#include <ipc/potentials/normal_potential.hpp>

namespace ipc {

/// @brief The normal adhesion potential.
class NormalAdhesionPotential : public NormalPotential {
    using Super = NormalPotential;

public:
    NormalAdhesionPotential(
        const double _dhat_p,
        const double _dhat_a,
        const double _Y,
        const double _eps_c);

    using Super::operator();
    using Super::gradient;
    using Super::hessian;

    /// @brief Compute the force magnitude for a collision.
    /// @param distance_squared The squared distance between elements.
    /// @param dmin The minimum distance offset to the barrier.
    /// @param barrier_stiffness The barrier stiffness.
    /// @return The force magnitude.
    double force_magnitude(
        const double distance_squared,
        const double dmin,
        const double barrier_stiffness) const override;

    /// @brief Compute the gradient of the force magnitude for a collision.
    /// @param distance_squared The squared distance between elements.
    /// @param distance_squared_gradient The gradient of the squared distance.
    /// @param dmin The minimum distance offset to the barrier.
    /// @param barrier_stiffness The stiffness of the barrier.
    /// @return The gradient of the force.
    VectorMax12d force_magnitude_gradient(
        const double distance_squared,
        const VectorMax12d& distance_squared_gradient,
        const double dmin,
        const double barrier_stiffness) const override;

    /// @brief The distance of largest adhesion force (\f$\hat{d}_p\f$) (\f$0 < \hat{d}_p < \hat{d}_a\f$).
    double dhat_p;
    /// @brief The adhesion activation distance (\f$\hat{d}_a\f$).
    double dhat_a;
    /// @brief The Young's modulus (\f$Y\f$).
    double Y;
    /// @brief The critical strain (\f$\varepsilon_c\f$).
    double eps_c;

protected:
    /// @brief Compute the barrier potential for a collision.
    /// @param distance_squared The distance (squared) between the two objects.
    /// @param dmin The minimum distance (unsquared) between the two objects.
    /// @return The barrier potential.
    double operator()(
        const double distance_squared, const double dmin = 0) const override;

    /// @brief Compute the gradient of the barrier potential for a collision.
    /// @param distance_squared The distance (squared) between the two objects.
    /// @param dmin The minimum distance (unsquared) between the two objects.
    /// @return The gradient of the barrier potential.
    double gradient(
        const double distance_squared, const double dmin = 0) const override;

    /// @brief Compute the hessian of the barrier potential for a collision.
    /// @param distance_squared The distance (squared) between the two objects.
    /// @param dmin The minimum distance (unsquared) between the two objects.
    /// @return The hessian of the barrier potential.
    double hessian(
        const double distance_squared, const double dmin = 0) const override;

private:
    std::array<double, 3>
    normal_adhesion_potential_args(const double dmin) const;
};

} // namespace ipc