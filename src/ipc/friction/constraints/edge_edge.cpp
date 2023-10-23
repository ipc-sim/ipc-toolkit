#include "edge_edge.hpp"

#include <ipc/friction/closest_point.hpp>
#include <ipc/friction/relative_velocity.hpp>
#include <ipc/friction/tangent_basis.hpp>
#include <ipc/distance/edge_edge.hpp>

namespace ipc {

EdgeEdgeFrictionConstraint::EdgeEdgeFrictionConstraint(
    const EdgeEdgeConstraint& constraint)
    : EdgeEdgeCandidate(constraint.edge0_id, constraint.edge1_id)
{
    this->weight = constraint.weight;
    this->weight_gradient = constraint.weight_gradient;
}

EdgeEdgeFrictionConstraint::EdgeEdgeFrictionConstraint(
    const EdgeEdgeConstraint& constraint,
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double dhat,
    const double barrier_stiffness)
    : EdgeEdgeFrictionConstraint(constraint)
{
    FrictionConstraint::init(
        vertices, edges, faces, dhat, barrier_stiffness, constraint.dmin);
}

// ============================================================================

MatrixMax<double, 3, 2> EdgeEdgeFrictionConstraint::compute_tangent_basis(
    const VectorMax12d& positions) const
{

    assert(positions.size() == ndof());
    return edge_edge_tangent_basis(
        positions.head(dim()), positions.segment(dim(), dim()),
        positions.segment(2 * dim(), dim()), positions.tail(dim()));
}

MatrixMax<double, 36, 2>
EdgeEdgeFrictionConstraint::compute_tangent_basis_jacobian(
    const VectorMax12d& positions) const
{
    assert(positions.size() == ndof());
    return edge_edge_tangent_basis_jacobian(
        positions.head(dim()), positions.segment(dim(), dim()),
        positions.segment(2 * dim(), dim()), positions.tail(dim()));
}

// ============================================================================

VectorMax2d EdgeEdgeFrictionConstraint::compute_closest_point(
    const VectorMax12d& positions) const
{
    assert(positions.size() == ndof());
    return edge_edge_closest_point(
        positions.head(dim()), positions.segment(dim(), dim()),
        positions.segment(2 * dim(), dim()), positions.tail(dim()));
}

MatrixMax<double, 2, 12>
EdgeEdgeFrictionConstraint::compute_closest_point_jacobian(
    const VectorMax12d& positions) const
{
    assert(positions.size() == ndof());
    return edge_edge_closest_point_jacobian(
        positions.head(dim()), positions.segment(dim(), dim()),
        positions.segment(2 * dim(), dim()), positions.tail(dim()));
}

// ============================================================================

VectorMax3d EdgeEdgeFrictionConstraint::relative_velocity(
    const VectorMax12d& velocities) const
{
    assert(velocities.size() == 12);
    return edge_edge_relative_velocity(
        velocities.head<3>(), velocities.segment<3>(dim()),
        velocities.segment<3>(2 * dim()), velocities.tail<3>(), closest_point);
}

MatrixMax<double, 3, 12> EdgeEdgeFrictionConstraint::relative_velocity_matrix(
    const VectorMax2d& _closest_point) const
{
    assert(_closest_point.size() == 2);
    return edge_edge_relative_velocity_matrix(dim(), _closest_point);
}

MatrixMax<double, 6, 12>
EdgeEdgeFrictionConstraint::relative_velocity_matrix_jacobian(
    const VectorMax2d& _closest_point) const
{
    assert(_closest_point.size() == 2);
    return edge_edge_relative_velocity_matrix_jacobian(dim(), _closest_point);
}

} // namespace ipc
