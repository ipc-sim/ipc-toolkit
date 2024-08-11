#include "barrier_potential.hpp"

#include <ipc/barrier/barrier.hpp>
#include <ipc/barrier/barrier_force_magnitude.hpp>

namespace ipc {

BarrierPotential::BarrierPotential(const double dhat) : NormalPotential()
{
    set_dhat(dhat);
}

BarrierPotential::BarrierPotential(
    const std::shared_ptr<Barrier> barrier, const double dhat)
    : NormalPotential()
{
    set_dhat(dhat);
    set_barrier(barrier);
}

double BarrierPotential::force_magnitude(
    const double distance_squared,
    const double dmin,
    const double barrier_stiffness) const
{
    return barrier_force_magnitude(
        distance_squared, barrier(), dhat(), barrier_stiffness, dmin);
}

VectorMax12d BarrierPotential::force_magnitude_gradient(
    const double distance_squared,
    const VectorMax12d& distance_squared_gradient,
    const double dmin,
    const double barrier_stiffness) const
{
    return barrier_force_magnitude_gradient(
        distance_squared, distance_squared_gradient, barrier(), dhat(),
        barrier_stiffness, dmin);
}

double BarrierPotential::operator()(
    const double distance_squared, const double dmin) const
{
    return barrier()(
        distance_squared - dmin * dmin, 2 * dmin * dhat() + dhat() * dhat());
}

double BarrierPotential::gradient(
    const double distance_squared, const double dmin) const
{
    return barrier().first_derivative(
        distance_squared - dmin * dmin, 2 * dmin * dhat() + dhat() * dhat());
}

double BarrierPotential::hessian(
    const double distance_squared, const double dmin) const
{
    return barrier().second_derivative(
        distance_squared - dmin * dmin, 2 * dmin * dhat() + dhat() * dhat());
}

} // namespace ipc