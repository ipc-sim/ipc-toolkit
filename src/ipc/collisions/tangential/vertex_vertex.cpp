#include "vertex_vertex.hpp"

#include <ipc/distance/point_point.hpp>
#include <ipc/tangent/closest_point.hpp>
#include <ipc/tangent/relative_velocity.hpp>
#include <ipc/tangent/tangent_basis.hpp>

namespace ipc {

VertexVertexTangentialCollision::VertexVertexTangentialCollision(
    const VertexVertexNormalCollision& collision)
    : VertexVertexCandidate(collision.vertex0_id, collision.vertex1_id)
{
    this->weight = collision.weight;
    this->weight_gradient = collision.weight_gradient;
    this->material_id1 = collision.material_id1;
    this->material_id2 = collision.material_id2;
}

VertexVertexTangentialCollision::VertexVertexTangentialCollision(
    const VertexVertexNormalCollision& collision,
    Eigen::ConstRef<VectorMax12d> positions,
    const NormalPotential& normal_potential,
    const double normal_stiffness)
    : VertexVertexTangentialCollision(collision)
{
    TangentialCollision::init(
        collision, positions, normal_potential, normal_stiffness);
}

// ============================================================================

MatrixMax<double, 3, 2> VertexVertexTangentialCollision::compute_tangent_basis(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    assert(positions.size() == ndof());
    return point_point_tangent_basis(
        positions.head(dim()), positions.tail(dim()));
}

MatrixMax<double, 36, 2>
VertexVertexTangentialCollision::compute_tangent_basis_jacobian(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    assert(positions.size() == ndof());
    return point_point_tangent_basis_jacobian(
        positions.head(dim()), positions.tail(dim()));
}

// ============================================================================

VectorMax2d VertexVertexTangentialCollision::compute_closest_point(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    return VectorMax2d();
}

MatrixMax<double, 2, 12>
VertexVertexTangentialCollision::compute_closest_point_jacobian(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    return MatrixMax<double, 2, 12>();
}

// ============================================================================

VectorMax3d VertexVertexTangentialCollision::relative_velocity(
    Eigen::ConstRef<VectorMax12d> velocities) const
{
    assert(velocities.size() == ndof());
    return point_point_relative_velocity(
        velocities.head(dim()), velocities.tail(dim()));
}

MatrixMax<double, 3, 12>
VertexVertexTangentialCollision::relative_velocity_matrix(
    Eigen::ConstRef<VectorMax2d> _closest_point) const
{
    assert(_closest_point.size() == 0);
    return point_point_relative_velocity_matrix(dim());
}

MatrixMax<double, 6, 12>
VertexVertexTangentialCollision::relative_velocity_matrix_jacobian(
    Eigen::ConstRef<VectorMax2d> _closest_point) const
{
    assert(_closest_point.size() == 0);
    return point_point_relative_velocity_matrix_jacobian(dim());
}

} // namespace ipc
