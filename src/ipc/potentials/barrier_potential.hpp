#pragma once

#include <ipc/potentials/distance_based_potential.hpp>
#include <ipc/barrier/barrier.hpp>

#include <memory>

namespace ipc {

/// @brief The barrier collision potential.
class BarrierPotential : public DistanceBasedPotential {
public:
    /// @brief Construct a barrier potential.
    /// @param dhat The activation distance of the barrier.
    explicit BarrierPotential(const double dhat);

    /// @brief Construct a barrier potential.
    /// @param dhat The activation distance of the barrier.
    BarrierPotential(const std::shared_ptr<Barrier> barrier, const double dhat);

    /// @brief Get the activation distance of the barrier.
    double dhat() const { return m_dhat; }

    /// @brief Set the activation distance of the barrier.
    /// @param dhat The activation distance of the barrier.
    void set_dhat(const double dhat)
    {
        assert(dhat > 0);
        m_dhat = dhat;
    }

    /// @brief Get the barrier function used to compute the potential.
    const Barrier& barrier() const
    {
        assert(m_barrier != nullptr);
        return *m_barrier;
    }

    /// @brief Set the barrier function used to compute the potential.
    /// @param barrier The barrier function used to compute the potential.
    void set_barrier(const std::shared_ptr<Barrier> barrier)
    {
        assert(barrier != nullptr);
        m_barrier = barrier;
    }

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

    /// @brief The activation distance of the barrier.
    double m_dhat;

    /// @brief The barrier function used to compute the potential.
    std::shared_ptr<Barrier> m_barrier = std::make_shared<ClampedLogBarrier>();
};

} // namespace ipc