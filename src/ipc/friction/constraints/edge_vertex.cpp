#include "edge_vertex.hpp"

#include <ipc/friction/closest_point.hpp>
#include <ipc/friction/relative_velocity.hpp>
#include <ipc/friction/tangent_basis.hpp>
#include <ipc/distance/point_edge.hpp>

namespace ipc {

EdgeVertexFrictionConstraint::EdgeVertexFrictionConstraint(
    long edge_index, long vertex_index)
    : EdgeVertexCandidate(edge_index, vertex_index)
{
}

EdgeVertexFrictionConstraint::EdgeVertexFrictionConstraint(
    const EdgeVertexCandidate& candidate)
    : EdgeVertexCandidate(candidate)
{
}

EdgeVertexFrictionConstraint::EdgeVertexFrictionConstraint(
    const EdgeVertexConstraint& constraint)
    : EdgeVertexCandidate(constraint.edge_index, constraint.vertex_index)
{
    this->weight = constraint.weight;
    this->weight_gradient = constraint.weight_gradient;
}

// ============================================================================

double
EdgeVertexFrictionConstraint::compute_distance(const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_edge_distance(
        x.head(dim()), x.segment(dim(), dim()), x.tail(dim()),
        PointEdgeDistanceType::P_E);
}

VectorMax12d EdgeVertexFrictionConstraint::compute_distance_gradient(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    VectorMax9d grad_d;
    point_edge_distance_gradient(
        x.head(dim()), x.segment(dim(), dim()), x.tail(dim()),
        PointEdgeDistanceType::P_E, grad_d);
    return grad_d;
}

// ============================================================================

MatrixMax<double, 3, 2>
EdgeVertexFrictionConstraint::compute_tangent_basis(const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_edge_tangent_basis(
        x.head(dim()), x.segment(dim(), dim()), x.tail(dim()));
}

MatrixMax<double, 36, 2>
EdgeVertexFrictionConstraint::compute_tangent_basis_jacobian(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_edge_tangent_basis_jacobian(
        x.head(dim()), x.segment(dim(), dim()), x.tail(dim()));
}

// ============================================================================

VectorMax2d
EdgeVertexFrictionConstraint::compute_closest_point(const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    VectorMax2d closest_point(1);
    closest_point[0] = point_edge_closest_point(
        x.head(dim()), x.segment(dim(), dim()), x.tail(dim()));
    return closest_point;
}

MatrixMax<double, 2, 12>
EdgeVertexFrictionConstraint::compute_closest_point_jacobian(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_edge_closest_point_jacobian(
               x.head(dim()), x.segment(dim(), dim()), x.tail(dim()))
        .transpose();
}

// ============================================================================

MatrixMax<double, 3, 12> EdgeVertexFrictionConstraint::relative_velocity_matrix(
    const VectorMax2d& closest_point) const
{
    assert(closest_point.size() == 1);
    return point_edge_relative_velocity_matrix(dim(), closest_point[0]);
}

MatrixMax<double, 6, 12>
EdgeVertexFrictionConstraint::relative_velocity_matrix_jacobian(
    const VectorMax2d& closest_point) const
{
    assert(closest_point.size() == 1);
    return point_edge_relative_velocity_matrix_jacobian(
        dim(), closest_point[0]);
}

} // namespace ipc
