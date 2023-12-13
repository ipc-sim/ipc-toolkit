#pragma once

#include <ipc/potentials/distance_based_potential.hpp>

namespace ipc {

class BarrierPotential : public DistanceBasedPotential {
public:
    /// @brief Construct a barrier potential.
    /// @param dhat The activation distance of the barrier.
    BarrierPotential(const double dhat);

    /// @brief Get the activation distance of the barrier.
    double dhat() const { return m_dhat; }

    /// @brief Set the activation distance of the barrier.
    /// @param dhat The activation distance of the barrier.
    void set_dhat(const double dhat)
    {
        assert(dhat > 0);
        m_dhat = dhat;
    }

protected:
    /// @brief Compute the barrier potential for a contact.
    /// @param distance_sqr The distance (squared) between the two objects.
    /// @param dmin The minimum distance (unsquared) between the two objects.
    /// @return The barrier potential.
    double distance_based_potential(
        const double distance_sqr, const double dmin = 0) const override;

    /// @brief Compute the gradient of the barrier potential for a contact.
    /// @param distance_sqr The distance (squared) between the two objects.
    /// @param dmin The minimum distance (unsquared) between the two objects.
    /// @return The gradient of the barrier potential.
    double distance_based_potential_gradient(
        const double distance_sqr, const double dmin = 0) const override;

    /// @brief Compute the hessian of the barrier potential for a contact.
    /// @param distance_sqr The distance (squared) between the two objects.
    /// @param dmin The minimum distance (unsquared) between the two objects.
    /// @return The hessian of the barrier potential.
    double distance_based_potential_hessian(
        const double distance_sqr, const double dmin = 0) const override;

    /// @brief The activation distance of the barrier.
    double m_dhat;
};

} // namespace ipc