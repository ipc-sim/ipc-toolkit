#include "vertex_vertex.hpp"

#include <ipc/friction/closest_point.hpp>
#include <ipc/friction/relative_velocity.hpp>
#include <ipc/friction/tangent_basis.hpp>
#include <ipc/distance/point_point.hpp>

namespace ipc {

VertexVertexFrictionConstraint::VertexVertexFrictionConstraint(
    const VertexVertexConstraint& constraint)
    : VertexVertexCandidate(constraint.vertex0_id, constraint.vertex1_id)
{
    this->weight = constraint.weight;
    this->weight_gradient = constraint.weight_gradient;
}

VertexVertexFrictionConstraint::VertexVertexFrictionConstraint(
    const VertexVertexConstraint& constraint,
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double dhat,
    const double barrier_stiffness)
    : VertexVertexFrictionConstraint(constraint)
{
    FrictionConstraint::init(
        vertices, edges, faces, dhat, barrier_stiffness, constraint.dmin);
}

// ============================================================================

MatrixMax<double, 3, 2> VertexVertexFrictionConstraint::compute_tangent_basis(
    const VectorMax12d& positions) const
{
    assert(positions.size() == ndof());
    return point_point_tangent_basis(
        positions.head(dim()), positions.tail(dim()));
}

MatrixMax<double, 36, 2>
VertexVertexFrictionConstraint::compute_tangent_basis_jacobian(
    const VectorMax12d& positions) const
{
    assert(positions.size() == ndof());
    return point_point_tangent_basis_jacobian(
        positions.head(dim()), positions.tail(dim()));
}

// ============================================================================

VectorMax2d VertexVertexFrictionConstraint::compute_closest_point(
    const VectorMax12d& positions) const
{
    return VectorMax2d();
}

MatrixMax<double, 2, 12>
VertexVertexFrictionConstraint::compute_closest_point_jacobian(
    const VectorMax12d& positions) const
{
    return MatrixMax<double, 2, 12>();
}

// ============================================================================

VectorMax3d VertexVertexFrictionConstraint::relative_velocity(
    const VectorMax12d& velocity) const
{
    assert(velocity.size() == ndof());
    return point_point_relative_velocity(
        velocity.head(dim()), velocity.tail(dim()));
}

MatrixMax<double, 3, 12>
VertexVertexFrictionConstraint::relative_velocity_matrix(
    const VectorMax2d& _closest_point) const
{
    assert(_closest_point.size() == 0);
    return point_point_relative_velocity_matrix(dim());
}

MatrixMax<double, 6, 12>
VertexVertexFrictionConstraint::relative_velocity_matrix_jacobian(
    const VectorMax2d& _closest_point) const
{
    assert(_closest_point.size() == 0);
    return point_point_relative_velocity_matrix_jacobian(dim());
}

} // namespace ipc
