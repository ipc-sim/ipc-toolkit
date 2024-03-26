#include "vertex_vertex.hpp"

#include <ipc/friction/closest_point.hpp>
#include <ipc/friction/relative_velocity.hpp>
#include <ipc/friction/tangent_basis.hpp>
#include <ipc/distance/point_point.hpp>

namespace ipc {

VertexVertexFrictionCollision::VertexVertexFrictionCollision(
    const VertexVertexCollision& collision)
    : VertexVertexCandidate(collision.vertex0_id, collision.vertex1_id)
{
    this->weight = collision.weight;
    this->weight_gradient = collision.weight_gradient;
}

VertexVertexFrictionCollision::VertexVertexFrictionCollision(
    const VertexVertexCollision& collision,
    const VectorMax12d& positions,
    const BarrierPotential& barrier_potential,
    const double barrier_stiffness)
    : VertexVertexFrictionCollision(collision)
{
    FrictionCollision::init(
        collision, positions, barrier_potential, barrier_stiffness);
}

// ============================================================================

MatrixMax<double, 3, 2> VertexVertexFrictionCollision::compute_tangent_basis(
    const VectorMax12d& positions) const
{
    assert(positions.size() == ndof());
    return point_point_tangent_basis(
        positions.head(dim()), positions.tail(dim()));
}

MatrixMax<double, 36, 2>
VertexVertexFrictionCollision::compute_tangent_basis_jacobian(
    const VectorMax12d& positions) const
{
    assert(positions.size() == ndof());
    return point_point_tangent_basis_jacobian(
        positions.head(dim()), positions.tail(dim()));
}

// ============================================================================

VectorMax2d VertexVertexFrictionCollision::compute_closest_point(
    const VectorMax12d& positions) const
{
    return VectorMax2d();
}

MatrixMax<double, 2, 12>
VertexVertexFrictionCollision::compute_closest_point_jacobian(
    const VectorMax12d& positions) const
{
    return MatrixMax<double, 2, 12>();
}

// ============================================================================

VectorMax3d VertexVertexFrictionCollision::relative_velocity(
    const VectorMax12d& velocities) const
{
    assert(velocities.size() == ndof());
    return point_point_relative_velocity(
        velocities.head(dim()), velocities.tail(dim()));
}

MatrixMax<double, 3, 12>
VertexVertexFrictionCollision::relative_velocity_matrix(
    const VectorMax2d& _closest_point) const
{
    assert(_closest_point.size() == 0);
    return point_point_relative_velocity_matrix(dim());
}

MatrixMax<double, 6, 12>
VertexVertexFrictionCollision::relative_velocity_matrix_jacobian(
    const VectorMax2d& _closest_point) const
{
    assert(_closest_point.size() == 0);
    return point_point_relative_velocity_matrix_jacobian(dim());
}

} // namespace ipc
