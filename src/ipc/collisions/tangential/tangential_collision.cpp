#include "tangential_collision.hpp"

#include <ipc/config.hpp>

namespace ipc {

void TangentialCollision::init(
    const NormalCollision& collision,
    Eigen::ConstRef<VectorMax12d> positions,
    const NormalPotential& normal_potential)
{
    init(
        collision, positions,
        normal_potential.force_magnitude(
            compute_distance(positions), collision.dmin));
}

void TangentialCollision::init(
    const NormalCollision& collision,
    Eigen::ConstRef<VectorMax12d> positions,
    const double normal_force)
{
    // do this to initialize dim()
    const int dim = collision.dim(positions.size());
    tangent_basis.resize(dim, dim - 1);

    closest_point = compute_closest_point(positions);
    tangent_basis = compute_tangent_basis(positions);
    normal_force_magnitude = normal_force;
}

} // namespace ipc
