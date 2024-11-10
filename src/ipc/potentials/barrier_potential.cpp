#include "barrier_potential.hpp"

#include <ipc/barrier/barrier.hpp>
#include <ipc/barrier/barrier_force_magnitude.hpp>

namespace ipc {

BarrierPotential::BarrierPotential(
    const double dhat, const bool use_physical_barrier)
    : BarrierPotential(
          std::make_shared<ClampedLogBarrier>(), dhat, use_physical_barrier)
{
}

BarrierPotential::BarrierPotential(
    const std::shared_ptr<Barrier> barrier,
    const double dhat,
    const bool use_physical_barrier)
    : NormalPotential()
{
    set_dhat(dhat);
    set_barrier(barrier);
    set_use_physical_barrier(use_physical_barrier);
}

double BarrierPotential::force_magnitude(
    const double distance_squared,
    const double dmin,
    const double barrier_stiffness) const
{
    double N = barrier_force_magnitude(
        distance_squared, barrier(), dhat(), barrier_stiffness, dmin);

    if (use_physical_barrier()) {
        N *= dhat() / barrier().units((2 * dmin + dhat()) * dhat());
    }

    return N;
}

VectorMax12d BarrierPotential::force_magnitude_gradient(
    const double distance_squared,
    const VectorMax12d& distance_squared_gradient,
    const double dmin,
    const double barrier_stiffness) const
{
    VectorMax12d grad_N = barrier_force_magnitude_gradient(
        distance_squared, distance_squared_gradient, barrier(), dhat(),
        barrier_stiffness, dmin);

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

    return b;
}

double BarrierPotential::gradient(
    const double distance_squared, const double dmin) const
{
    double db = barrier().first_derivative(
        distance_squared - dmin * dmin, (2 * dmin + dhat()) * dhat());

    if (use_physical_barrier()) {
        db *= dhat() / barrier().units((2 * dmin + dhat()) * dhat());
    }

    return db;
}

double BarrierPotential::hessian(
    const double distance_squared, const double dmin) const
{
    double d2b = barrier().second_derivative(
        distance_squared - dmin * dmin, (2 * dmin + dhat()) * dhat());

    if (use_physical_barrier()) {
        d2b *= dhat() / barrier().units((2 * dmin + dhat()) * dhat());
    }

    return d2b;
}

} // namespace ipc