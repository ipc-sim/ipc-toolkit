#pragma once

#include <ipc/barrier/barrier.hpp>
#include <ipc/potentials/normal_potential.hpp>

#include <memory>

namespace ipc {

/// @brief The barrier collision potential.
class BarrierPotential : public NormalPotential {
    using Super = NormalPotential;

    // Delete the constructor that doesn't take stiffness to avoid
    // accidentally using the default stiffness value of 0.0, which
    // is likely not what we want.
    BarrierPotential(const double, const bool) = delete;

public:
    /// @brief Construct a barrier potential.
    /// @param dhat The activation distance of the barrier.
    /// @param stiffness The stiffness of the barrier.
    /// @param use_physical_barrier Whether to use the physical barrier.
    BarrierPotential(
        const double dhat,
        const double stiffness,
        const bool use_physical_barrier = false);

    /// @brief Construct a barrier potential.
    /// @param barrier The barrier function.
    /// @param dhat The activation distance of the barrier.
    /// @param stiffness The stiffness of the barrier.
    /// @param use_physical_barrier Whether to use the physical barrier.
    BarrierPotential(
        std::shared_ptr<Barrier> barrier,
        const double dhat,
        const double stiffness,
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

    /// @brief Get the stiffness of the barrier.
    double stiffness() const { return m_stiffness; }

    /// @brief Set the stiffness of the barrier.
    /// @param stiffness The stiffness of the barrier.
    void set_stiffness(const double stiffness)
    {
        assert(stiffness > 0);
        m_stiffness = stiffness;
    }

    /// @brief Get the barrier function used to compute the potential.
    const Barrier& barrier() const
    {
        assert(m_barrier != nullptr);
        return *m_barrier;
    }

    /// @brief Set the barrier function used to compute the potential.
    /// @param barrier The barrier function used to compute the potential.
    void set_barrier(const std::shared_ptr<Barrier>& barrier)
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
    /// @return The force magnitude.
    double force_magnitude(
        const double distance_squared, const double dmin) const override;

    /// @brief Compute the gradient of the force magnitude for a collision.
    /// @param distance_squared The squared distance between elements.
    /// @param distance_squared_gradient The gradient of the squared distance.
    /// @param dmin The minimum distance offset to the barrier.
    /// @return The gradient of the force.
    VectorMax12d force_magnitude_gradient(
        const double distance_squared,
        Eigen::ConstRef<VectorMax12d> distance_squared_gradient,
        const double dmin) const override;

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

    /// @brief The stiffness of the barrier.
    double m_stiffness;

    /// @brief Whether to use the physical barrier.
    /// @note When using the convergent formulation we want the barrier to
    ///       have units of Pa⋅m, so κ gets units of Pa and the barrier function
    ///       should have units of m. See notebooks/physical_barrier.ipynb for
    ///       more details.
    bool m_use_physical_barrier = false;
};

} // namespace ipc