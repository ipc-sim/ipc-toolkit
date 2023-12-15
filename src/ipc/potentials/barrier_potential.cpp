#include "barrier_potential.hpp"

#include <ipc/barrier/barrier.hpp>

namespace ipc {

BarrierPotential::BarrierPotential(const double dhat) : DistanceBasedPotential()
{
    set_dhat(dhat);
}

double BarrierPotential::distance_based_potential(
    const double distance_sqr, const double dmin) const
{
    assert(barrier != nullptr);
    return (*barrier)(
        distance_sqr - dmin * dmin, 2 * dmin * dhat() + dhat() * dhat());
}

double BarrierPotential::distance_based_potential_gradient(
    const double distance_sqr, const double dmin) const
{
    assert(barrier != nullptr);
    return barrier->first_derivative(
        distance_sqr - dmin * dmin, 2 * dmin * dhat() + dhat() * dhat());
}

double BarrierPotential::distance_based_potential_hessian(
    const double distance_sqr, const double dmin) const
{
    assert(barrier != nullptr);
    return barrier->second_derivative(
        distance_sqr - dmin * dmin, 2 * dmin * dhat() + dhat() * dhat());
}

} // namespace ipc