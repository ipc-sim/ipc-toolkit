#include "tangential_collision.hpp"

#include <ipc/config.hpp>

#include <stdexcept> // std::out_of_range

namespace ipc {

void TangentialCollision::init(
    const NormalCollision& collision,
    const VectorMax12d& positions,
    const NormalPotential& normal_potential,
    const double normal_stiffness)
{
    // do this to initialize dim()
    const int dim = collision.dim(positions.size());
    tangent_basis.resize(dim, dim - 1);

    closest_point = compute_closest_point(positions);
    tangent_basis = compute_tangent_basis(positions);
    normal_force_magnitude = normal_potential.force_magnitude(
        compute_distance(positions), collision.dmin, normal_stiffness);
}

} // namespace ipc
