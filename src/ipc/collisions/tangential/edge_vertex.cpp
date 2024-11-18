#include "edge_vertex.hpp"

#include <ipc/distance/point_edge.hpp>
#include <ipc/tangent/closest_point.hpp>
#include <ipc/tangent/relative_velocity.hpp>
#include <ipc/tangent/tangent_basis.hpp>

namespace ipc {

EdgeVertexTangentialCollision::EdgeVertexTangentialCollision(
    const EdgeVertexNormalCollision& collision)
    : EdgeVertexCandidate(collision.edge_id, collision.vertex_id)
{
    this->weight = collision.weight;
    this->weight_gradient = collision.weight_gradient;
}

EdgeVertexTangentialCollision::EdgeVertexTangentialCollision(
    const EdgeVertexNormalCollision& collision,
    const VectorMax12d& positions,
    const NormalPotential& normal_potential,
    const double normal_stiffness)
    : EdgeVertexTangentialCollision(collision)
{
    TangentialCollision::init(
        collision, positions, normal_potential, normal_stiffness);
}

// ============================================================================

MatrixMax<double, 3, 2> EdgeVertexTangentialCollision::compute_tangent_basis(
    const VectorMax12d& positions) const
{
    assert(positions.size() == ndof());
    return point_edge_tangent_basis(
        positions.head(dim()), positions.segment(dim(), dim()),
        positions.tail(dim()));
}

MatrixMax<double, 36, 2>
EdgeVertexTangentialCollision::compute_tangent_basis_jacobian(
    const VectorMax12d& positions) const
{
    assert(positions.size() == ndof());
    return point_edge_tangent_basis_jacobian(
        positions.head(dim()), positions.segment(dim(), dim()),
        positions.tail(dim()));
}

// ============================================================================

VectorMax2d EdgeVertexTangentialCollision::compute_closest_point(
    const VectorMax12d& positions) const
{
    assert(positions.size() == ndof());
    VectorMax2d alpha(1);
    alpha[0] = point_edge_closest_point(
        positions.head(dim()), positions.segment(dim(), dim()),
        positions.tail(dim()));
    return alpha;
}

MatrixMax<double, 2, 12>
EdgeVertexTangentialCollision::compute_closest_point_jacobian(
    const VectorMax12d& positions) const
{
    assert(positions.size() == ndof());
    return point_edge_closest_point_jacobian(
               positions.head(dim()), positions.segment(dim(), dim()),
               positions.tail(dim()))
        .transpose();
}

// ============================================================================

VectorMax3d EdgeVertexTangentialCollision::relative_velocity(
    const VectorMax12d& velocities) const
{
    assert(velocities.size() == ndof());
    return point_edge_relative_velocity(
        velocities.head(dim()), velocities.segment(dim(), dim()),
        velocities.tail(dim()), closest_point[0]);
}

MatrixMax<double, 3, 12>
EdgeVertexTangentialCollision::relative_velocity_matrix(
    const VectorMax2d& _closest_point) const
{
    assert(_closest_point.size() == 1);
    return point_edge_relative_velocity_matrix(dim(), _closest_point[0]);
}

MatrixMax<double, 6, 12>
EdgeVertexTangentialCollision::relative_velocity_matrix_jacobian(
    const VectorMax2d& _closest_point) const
{
    assert(_closest_point.size() == 1);
    return point_edge_relative_velocity_matrix_jacobian(
        dim(), _closest_point[0]);
}

} // namespace ipc
