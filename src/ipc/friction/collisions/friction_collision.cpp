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

    // Initialize dimension and tangent_basis based on the collision dimension
    const int dim = collision.dim(positions.size());
    tangent_basis.resize(dim, dim - 1);

    // Set the static and kinetic friction coefficients
    static_mu = 0.5;
    kinetic_mu = 0.5;

    closest_point = compute_closest_point(positions);
    tangent_basis = compute_tangent_basis(positions);

    // Compute the normal force magnitude based on barrier potential and stiffness
    normal_force_magnitude = this->compute_normal_force_magnitude(
        positions, barrier_potential, barrier_stiffness, collision.dmin);

}

void FrictionCollision::init(
    const Collision& collision,
    const VectorMax12d& positions,
    const BarrierPotential& barrier_potential,
    const double barrier_stiffness,
    double static_mu,
    double kinetic_mu)
{
    // Initialize dimension and tangent_basis based on the collision dimension
    const int dim = collision.dim(positions.size());
    tangent_basis.resize(dim, dim - 1);

    // Set the static and kinetic friction coefficients
    this->static_mu = static_mu;
    this->kinetic_mu = kinetic_mu;

    closest_point = compute_closest_point(positions);
    tangent_basis = compute_tangent_basis(positions);

    // Compute the normal force magnitude based on barrier potential and stiffness
    normal_force_magnitude = this->compute_normal_force_magnitude(
        positions, barrier_potential, barrier_stiffness, collision.dmin);
}

double FrictionCollision::compute_normal_force_magnitude(
    const VectorMax12d& positions,
    const BarrierPotential& barrier_potential,
    const double barrier_stiffness,
    const double dmin) const
{
    const double dhat = barrier_potential.dhat();
    double N = ipc::compute_normal_force_magnitude(
        compute_distance(positions), barrier_potential.barrier(), dhat,
        barrier_stiffness, dmin);

    // Adjust normal force magnitude for physical barrier if necessary
    if (barrier_potential.use_physical_barrier()) {
        N *= dhat / barrier_potential.barrier().units((2 * dmin + dhat) * dhat);
    }

    return N;
}

VectorMax12d FrictionCollision::compute_normal_force_magnitude_gradient(
    const VectorMax12d& positions,
    const BarrierPotential& barrier_potential,
    const double barrier_stiffness,
    const double dmin) const
{
    const double dhat = barrier_potential.dhat();
    VectorMax12d grad_N = ipc::compute_normal_force_magnitude_gradient(
        compute_distance(positions), compute_distance_gradient(positions),
        barrier_potential.barrier(), dhat, barrier_stiffness, dmin);

    // Adjust gradient of the normal force magnitude if physical barrier is used
    if (barrier_potential.use_physical_barrier()) {
        grad_N *= dhat / barrier_potential.barrier().units((2 * dmin + dhat) * dhat);
    }

    return grad_N;
}


} // namespace ipc
