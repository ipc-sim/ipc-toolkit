#include "barrier_potential.hpp"

#include <ipc/barrier/barrier.hpp>

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
    : DistanceBasedPotential()
{
    set_dhat(dhat);
    set_barrier(barrier);
    set_use_physical_barrier(use_physical_barrier);
}

double BarrierPotential::distance_based_potential(
    const double distance_sqr, const double dmin) const
{
    double b =
        barrier()(distance_sqr - dmin * dmin, (2 * dmin + dhat()) * dhat());

    if (use_physical_barrier()) {
        b *= dhat() / barrier().units((2 * dmin + dhat()) * dhat());
    }

    return b;
}

double BarrierPotential::distance_based_potential_gradient(
    const double distance_sqr, const double dmin) const
{
    double db = barrier().first_derivative(
        distance_sqr - dmin * dmin, (2 * dmin + dhat()) * dhat());

    if (use_physical_barrier()) {
        db *= dhat() / barrier().units((2 * dmin + dhat()) * dhat());
    }

    return db;
}

double BarrierPotential::distance_based_potential_hessian(
    const double distance_sqr, const double dmin) const
{
    double d2b = barrier().second_derivative(
        distance_sqr - dmin * dmin, (2 * dmin + dhat()) * dhat());

    if (use_physical_barrier()) {
        d2b *= dhat() / barrier().units((2 * dmin + dhat()) * dhat());
    }

    return d2b;
}

} // namespace ipc