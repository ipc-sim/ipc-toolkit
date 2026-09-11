#include "barrier_potential.hpp"

#include <ipc/barrier/barrier.hpp>
#include <ipc/barrier/barrier_force_magnitude.hpp>

#include <cmath>
#include <stdexcept>

namespace ipc {

BarrierPotential::BarrierPotential(
    const double dhat, const double stiffness, const bool use_physical_barrier)
    : BarrierPotential(
          std::make_shared<ClampedLogBarrier<>>(),
          dhat,
          stiffness,
          use_physical_barrier)
{
}

BarrierPotential::BarrierPotential(
    std::shared_ptr<Barrier> barrier,
    const double dhat,
    const double stiffness,
    const bool use_physical_barrier,
    const bool use_squared_distance)
    : m_barrier(std::move(barrier))
    , m_dhat(dhat)
    , m_stiffness(stiffness)
    , m_use_physical_barrier(use_physical_barrier)
    , m_use_squared_distance(use_squared_distance)
{
    assert(dhat > 0);
    assert(stiffness > 0);
    assert(m_barrier != nullptr);
}

double BarrierPotential::force_magnitude(
    const double distance_squared, const double dmin) const
{
    if (!m_use_squared_distance) {
        throw std::runtime_error(
            "BarrierPotential: force_magnitude not implemented for "
            "use_squared_distance=false");
    }

    double N = barrier_force_magnitude(
        distance_squared, barrier(), dhat(), stiffness(), dmin);

    if (use_physical_barrier()) {
        N *= dhat() / barrier().units((2 * dmin + dhat()) * dhat());
    }

    return N;
}

VectorMax12d BarrierPotential::force_magnitude_gradient(
    const double distance_squared,
    Eigen::ConstRef<VectorMax12d> distance_squared_gradient,
    const double dmin) const
{
    if (!m_use_squared_distance) {
        throw std::runtime_error(
            "BarrierPotential: force_magnitude_gradient not implemented for "
            "use_squared_distance=false");
    }

    VectorMax12d grad_N = barrier_force_magnitude_gradient(
        distance_squared, distance_squared_gradient, barrier(), dhat(),
        stiffness(), dmin);

    if (use_physical_barrier()) {
        grad_N *= dhat() / barrier().units((2 * dmin + dhat()) * dhat());
    }

    return grad_N;
}

double BarrierPotential::operator()(
    const double distance_squared, const double dmin) const
{
    if (!m_use_squared_distance) {
        const double d = std::sqrt(distance_squared);
        double b = barrier()(d - dmin, dhat());
        if (use_physical_barrier()) {
            b *= dhat() / barrier().units(dhat());
        }
        return b;
    }

    double b =
        barrier()(distance_squared - dmin * dmin, (2 * dmin + dhat()) * dhat());

    if (use_physical_barrier()) {
        b *= dhat() / barrier().units((2 * dmin + dhat()) * dhat());
    }

    return stiffness() * b;
}

double BarrierPotential::gradient(
    const double distance_squared, const double dmin) const
{
    if (!m_use_squared_distance) {
        const double d = std::sqrt(distance_squared);
        double db = barrier().first_derivative(d - dmin, dhat()) / (2.0 * d);
        if (use_physical_barrier()) {
            db *= dhat() / barrier().units(dhat());
        }
        return db;
    }

    double db = barrier().first_derivative(
        distance_squared - dmin * dmin, (2 * dmin + dhat()) * dhat());

    if (use_physical_barrier()) {
        db *= dhat() / barrier().units((2 * dmin + dhat()) * dhat());
    }

    return stiffness() * db;
}

double BarrierPotential::hessian(
    const double distance_squared, const double dmin) const
{
    if (!m_use_squared_distance) {
        const double d = std::sqrt(distance_squared);
        const double b1 = barrier().first_derivative(d - dmin, dhat());
        const double b2 = barrier().second_derivative(d - dmin, dhat());
        double d2b =
            b2 / (4.0 * distance_squared) - b1 / (4.0 * d * distance_squared);
        if (use_physical_barrier()) {
            d2b *= dhat() / barrier().units(dhat());
        }
        return d2b;
    }

    double d2b = barrier().second_derivative(
        distance_squared - dmin * dmin, (2 * dmin + dhat()) * dhat());

    if (use_physical_barrier()) {
        d2b *= dhat() / barrier().units((2 * dmin + dhat()) * dhat());
    }

    return stiffness() * d2b;
}

} // namespace ipc