#include "tangential_collision.hpp"

#include <ipc/config.hpp>

#include <stdexcept> // std::out_of_range

namespace ipc {

void TangentialCollision::init(
    const NormalCollision& collision,
    Eigen::ConstRef<VectorMax12d> positions,
    const NormalPotential& normal_potential,
    const double normal_stiffness)
{
    this->weight = collision.weight;
    this->weight_gradient = collision.weight_gradient;
    this->material_id1 = collision.material_id1;
    this->material_id2 = collision.material_id2;

    this->tangent_basis = compute_tangent_basis(positions);
    this->closest_point = compute_closest_point(positions);
    double distance = collision.compute_distance(positions);

    // We should compute a more accurate distance mollifier here.
    this->normal_force_magnitude = normal_potential.force_magnitude(
        distance, collision.dmin, normal_stiffness);
}

MatrixMax<double, 3, 12> TangentialCollision::relative_velocity_matrix() const
{
    return relative_velocity_matrix(closest_point);
}

int TangentialCollision::dim() const
{
    return 3;
}

int TangentialCollision::ndof() const
{
    return num_vertices() * dim();
}

int TangentialCollision::num_vertices() const
{
    return 0;
}

} // namespace ipc
