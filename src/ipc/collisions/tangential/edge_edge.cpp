#include "edge_edge.hpp"

#include <ipc/distance/edge_edge.hpp>
#include <ipc/tangent/closest_point.hpp>
#include <ipc/tangent/relative_velocity.hpp>
#include <ipc/tangent/tangent_basis.hpp>

namespace ipc {

EdgeEdgeTangentialCollision::EdgeEdgeTangentialCollision(
    const EdgeEdgeNormalCollision& collision)
    : EdgeEdgeCandidate(collision.edge0_id, collision.edge1_id)
{
    this->weight = collision.weight;
    this->weight_gradient = collision.weight_gradient;
}

EdgeEdgeTangentialCollision::EdgeEdgeTangentialCollision(
    const EdgeEdgeNormalCollision& collision,
    const VectorMax12d& positions,
    const NormalPotential& normal_potential,
    const double normal_stiffness)
    : EdgeEdgeTangentialCollision(collision)
{
    TangentialCollision::init(
        collision, positions, normal_potential, normal_stiffness);
}

// ============================================================================

MatrixMax<double, 3, 2> EdgeEdgeTangentialCollision::compute_tangent_basis(
    const VectorMax12d& positions) const
{

    assert(positions.size() == ndof());
    return edge_edge_tangent_basis(
        positions.head(dim()), positions.segment(dim(), dim()),
        positions.segment(2 * dim(), dim()), positions.tail(dim()));
}

MatrixMax<double, 36, 2>
EdgeEdgeTangentialCollision::compute_tangent_basis_jacobian(
    const VectorMax12d& positions) const
{
    assert(positions.size() == ndof());
    return edge_edge_tangent_basis_jacobian(
        positions.head(dim()), positions.segment(dim(), dim()),
        positions.segment(2 * dim(), dim()), positions.tail(dim()));
}

// ============================================================================

VectorMax2d EdgeEdgeTangentialCollision::compute_closest_point(
    const VectorMax12d& positions) const
{
    assert(positions.size() == ndof());
    return edge_edge_closest_point(
        positions.head(dim()), positions.segment(dim(), dim()),
        positions.segment(2 * dim(), dim()), positions.tail(dim()));
}

MatrixMax<double, 2, 12>
EdgeEdgeTangentialCollision::compute_closest_point_jacobian(
    const VectorMax12d& positions) const
{
    assert(positions.size() == ndof());
    return edge_edge_closest_point_jacobian(
        positions.head(dim()), positions.segment(dim(), dim()),
        positions.segment(2 * dim(), dim()), positions.tail(dim()));
}

// ============================================================================

VectorMax3d EdgeEdgeTangentialCollision::relative_velocity(
    const VectorMax12d& velocities) const
{
    assert(velocities.size() == 12);
    return edge_edge_relative_velocity(
        velocities.head<3>(), velocities.segment<3>(dim()),
        velocities.segment<3>(2 * dim()), velocities.tail<3>(), closest_point);
}

MatrixMax<double, 3, 12> EdgeEdgeTangentialCollision::relative_velocity_matrix(
    const VectorMax2d& _closest_point) const
{
    assert(_closest_point.size() == 2);
    return edge_edge_relative_velocity_matrix(dim(), _closest_point);
}

MatrixMax<double, 6, 12>
EdgeEdgeTangentialCollision::relative_velocity_matrix_jacobian(
    const VectorMax2d& _closest_point) const
{
    assert(_closest_point.size() == 2);
    return edge_edge_relative_velocity_matrix_jacobian(dim(), _closest_point);
}

} // namespace ipc
