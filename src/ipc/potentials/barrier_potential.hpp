#pragma once

#include <ipc/barrier/barrier.hpp>
#include <ipc/potentials/normal_potential.hpp>

#include <memory>

namespace ipc {

/// @brief The barrier collision potential.
class BarrierPotential : public NormalPotential {
    using Super = NormalPotential;

public:
    /// @brief Construct a barrier potential.
    /// @param dhat The activation distance of the barrier.
    explicit BarrierPotential(
        const double dhat, const bool use_physical_barrier = false);

    /// @brief Construct a barrier potential.
    /// @param barrier The barrier function.
    /// @param dhat The activation distance of the barrier.
    BarrierPotential(
        const std::shared_ptr<Barrier> barrier,
        const double dhat,
        const bool use_physical_barrier = false);

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

    /// @brief Get whether to use the physical barrier.
    /// @note When using the convergent formulation we want the barrier to
    ///       have units of Pa⋅m, so κ gets units of Pa and the barrier function
    ///       should have units of m. See notebooks/physical_barrier.ipynb for
    ///       more details.
    bool use_physical_barrier() const { return m_use_physical_barrier; }

    /// @brief Set use physical barrier flag.
    /// @param use_physical_barrier Whether to use the physical barrier.
    /// @note When using the convergent formulation we want the barrier to
    ///       have units of Pa⋅m, so κ gets units of Pa and the barrier function
    ///       should have units of m. See notebooks/physical_barrier.ipynb for
    ///       more details.
    void set_use_physical_barrier(bool use_physical_barrier)
    {
        m_use_physical_barrier = use_physical_barrier;
    }

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

    /// @brief The barrier function used to compute the potential.
    std::shared_ptr<Barrier> m_barrier = std::make_shared<ClampedLogBarrier>();

    /// @brief The activation distance of the barrier.
    double m_dhat;

    /// @brief Whether to use the physical barrier.
    /// @note When using the convergent formulation we want the barrier to
    ///       have units of Pa⋅m, so κ gets units of Pa and the barrier function
    ///       should have units of m. See notebooks/physical_barrier.ipynb for
    ///       more details.
    bool m_use_physical_barrier = false;
};

} // namespace ipc