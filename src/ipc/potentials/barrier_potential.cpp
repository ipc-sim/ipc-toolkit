#include "barrier_potential.hpp"

#include <ipc/barrier/barrier.hpp>
#include <ipc/barrier/barrier_force_magnitude.hpp>

namespace ipc {

BarrierPotential::BarrierPotential(
    const double dhat, const double stiffness, const bool use_physical_barrier)
    : BarrierPotential(
          std::make_shared<ClampedLogBarrier>(),
          dhat,
          stiffness,
          use_physical_barrier)
{
}

BarrierPotential::BarrierPotential(
    std::shared_ptr<Barrier> barrier,
    const double dhat,
    const double stiffness,
    const bool use_physical_barrier)
    : m_barrier(std::move(barrier))
    , m_dhat(dhat)
    , m_stiffness(stiffness)
    , m_use_physical_barrier(use_physical_barrier)
{
    assert(dhat > 0);
    assert(stiffness > 0);
    assert(m_barrier != nullptr);
}

double BarrierPotential::force_magnitude(
    const double distance_squared, const double dmin) const
{
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
    double d2b = barrier().second_derivative(
        distance_squared - dmin * dmin, (2 * dmin + dhat()) * dhat());

    if (use_physical_barrier()) {
        d2b *= dhat() / barrier().units((2 * dmin + dhat()) * dhat());
    }

    return stiffness() * d2b;
}

} // namespace ipc