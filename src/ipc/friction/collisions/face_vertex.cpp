#include "face_vertex.hpp"

#include <ipc/distance/point_triangle.hpp>
#include <ipc/friction/closest_point.hpp>
#include <ipc/friction/relative_velocity.hpp>
#include <ipc/friction/tangent_basis.hpp>

namespace ipc {

FaceVertexFrictionCollision::FaceVertexFrictionCollision(
    const FaceVertexCollision& collision)
    : FaceVertexCandidate(collision.face_id, collision.vertex_id)
{
    this->weight = collision.weight;
    this->weight_gradient = collision.weight_gradient;
}

FaceVertexFrictionCollision::FaceVertexFrictionCollision(
    const FaceVertexCollision& collision,
    const VectorMax12d& positions,
    const BarrierPotential& barrier_potential,
    const double barrier_stiffness)
    : FaceVertexFrictionCollision(collision)
{
    FrictionCollision::init(
        collision, positions, barrier_potential, barrier_stiffness);
}

// ============================================================================

MatrixMax<double, 3, 2> FaceVertexFrictionCollision::compute_tangent_basis(
    const VectorMax12d& positions) const
{
    assert(positions.size() == ndof());
    return point_triangle_tangent_basis(
        positions.head(dim()), positions.segment(dim(), dim()),
        positions.segment(2 * dim(), dim()), positions.tail(dim()));
}

MatrixMax<double, 36, 2>
FaceVertexFrictionCollision::compute_tangent_basis_jacobian(
    const VectorMax12d& positions) const
{
    assert(positions.size() == ndof());
    return point_triangle_tangent_basis_jacobian(
        positions.head(dim()), positions.segment(dim(), dim()),
        positions.segment(2 * dim(), dim()), positions.tail(dim()));
}

// ============================================================================

VectorMax2d FaceVertexFrictionCollision::compute_closest_point(
    const VectorMax12d& positions) const
{
    assert(positions.size() == ndof());
    return point_triangle_closest_point(
        positions.head(dim()), positions.segment(dim(), dim()),
        positions.segment(2 * dim(), dim()), positions.tail(dim()));
}

MatrixMax<double, 2, 12>
FaceVertexFrictionCollision::compute_closest_point_jacobian(
    const VectorMax12d& positions) const
{
    assert(positions.size() == ndof());
    return point_triangle_closest_point_jacobian(
        positions.head(dim()), positions.segment(dim(), dim()),
        positions.segment(2 * dim(), dim()), positions.tail(dim()));
}

// ============================================================================

VectorMax3d FaceVertexFrictionCollision::relative_velocity(
    const VectorMax12d& velocities) const
{
    assert(velocities.size() == 12);
    return point_triangle_relative_velocity(
        velocities.head<3>(), velocities.segment<3>(dim()),
        velocities.segment<3>(2 * dim()), velocities.tail<3>(), closest_point);
}

MatrixMax<double, 3, 12> FaceVertexFrictionCollision::relative_velocity_matrix(
    const VectorMax2d& _closest_point) const
{
    assert(_closest_point.size() == 2);
    return point_triangle_relative_velocity_matrix(dim(), _closest_point);
}

MatrixMax<double, 6, 12>
FaceVertexFrictionCollision::relative_velocity_matrix_jacobian(
    const VectorMax2d& _closest_point) const
{
    assert(_closest_point.size() == 2);
    return point_triangle_relative_velocity_matrix_jacobian(
        dim(), _closest_point);
}

} // namespace ipc
