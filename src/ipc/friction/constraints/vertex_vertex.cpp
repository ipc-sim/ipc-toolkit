#include "vertex_vertex.hpp"

#include <ipc/friction/closest_point.hpp>
#include <ipc/friction/relative_velocity.hpp>
#include <ipc/friction/tangent_basis.hpp>
#include <ipc/distance/point_point.hpp>

namespace ipc {

VertexVertexFrictionConstraint::VertexVertexFrictionConstraint(
    long vertex0_id, long vertex1_id)
    : VertexVertexCandidate(vertex0_id, vertex1_id)
{
}

VertexVertexFrictionConstraint::VertexVertexFrictionConstraint(
    const VertexVertexCandidate& candidate)
    : VertexVertexCandidate(candidate)
{
}

VertexVertexFrictionConstraint::VertexVertexFrictionConstraint(
    const VertexVertexConstraint& constraint)
    : VertexVertexCandidate(constraint.vertex0_id, constraint.vertex1_id)
{
    this->weight = constraint.weight;
    this->weight_gradient = constraint.weight_gradient;
}

// ============================================================================

double
VertexVertexFrictionConstraint::compute_distance(const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_point_distance(x.head(dim()), x.tail(dim()));
}

VectorMax12d VertexVertexFrictionConstraint::compute_distance_gradient(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    VectorMax6d grad_d;
    point_point_distance_gradient(x.head(dim()), x.tail(dim()), grad_d);
    return grad_d;
}

// ============================================================================

MatrixMax<double, 3, 2> VertexVertexFrictionConstraint::compute_tangent_basis(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_point_tangent_basis(x.head(dim()), x.tail(dim()));
}

MatrixMax<double, 36, 2>
VertexVertexFrictionConstraint::compute_tangent_basis_jacobian(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_point_tangent_basis_jacobian(x.head(dim()), x.tail(dim()));
}

// ============================================================================

VectorMax2d VertexVertexFrictionConstraint::compute_closest_point(
    const VectorMax12d& x) const
{
    return VectorMax2d();
}

MatrixMax<double, 2, 12>
VertexVertexFrictionConstraint::compute_closest_point_jacobian(
    const VectorMax12d& x) const
{
    return MatrixMax<double, 2, 12>();
}

// ============================================================================

MatrixMax<double, 3, 12>
VertexVertexFrictionConstraint::relative_velocity_matrix(
    const VectorMax2d& closest_point) const
{
    assert(closest_point.size() == 0);
    return point_point_relative_velocity_matrix(dim());
}

MatrixMax<double, 6, 12>
VertexVertexFrictionConstraint::relative_velocity_matrix_jacobian(
    const VectorMax2d& closest_point) const
{
    assert(closest_point.size() == 0);
    return point_point_relative_velocity_matrix_jacobian(dim());
}

} // namespace ipc
