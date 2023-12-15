#include "friction_collision.hpp"

#include <ipc/friction/normal_force_magnitude.hpp>

#include <ipc/config.hpp>

#include <stdexcept> // std::out_of_range

namespace ipc {

void FrictionCollision::init(
    const Collision& collision,
    const VectorMax12d& positions,
    const BarrierPotential& barrier_potential,
    const double barrier_stiffness)
{
    // do this to initialize dim()
    const int dim = collision.dim(positions.size());
    tangent_basis.resize(dim, dim - 1);

    closest_point = compute_closest_point(positions);
    tangent_basis = compute_tangent_basis(positions);
    normal_force_magnitude = this->compute_normal_force_magnitude(
        positions, barrier_potential, barrier_stiffness, collision.dmin);
}

double FrictionCollision::compute_normal_force_magnitude(
    const VectorMax12d& positions,
    const BarrierPotential& barrier_potential,
    const double barrier_stiffness,
    const double dmin) const
{
    return ipc::compute_normal_force_magnitude(
        compute_distance(positions), barrier_potential.barrier(),
        barrier_potential.dhat(), barrier_stiffness, dmin);
}

VectorMax12d FrictionCollision::compute_normal_force_magnitude_gradient(
    const VectorMax12d& positions,
    const BarrierPotential& barrier_potential,
    const double barrier_stiffness,
    const double dmin) const
{
    return ipc::compute_normal_force_magnitude_gradient(
        compute_distance(positions), compute_distance_gradient(positions),
        barrier_potential.barrier(), barrier_potential.dhat(),
        barrier_stiffness, dmin);
}

} // namespace ipc
