#include "edge_edge.hpp"

#include <ipc/friction/closest_point.hpp>
#include <ipc/friction/relative_velocity.hpp>
#include <ipc/friction/tangent_basis.hpp>
#include <ipc/distance/edge_edge.hpp>

namespace ipc {

EdgeEdgeFrictionConstraint::EdgeEdgeFrictionConstraint(
    long edge0_id, long edge1_id)
    : EdgeEdgeCandidate(edge0_id, edge1_id)
{
}

EdgeEdgeFrictionConstraint::EdgeEdgeFrictionConstraint(
    const EdgeEdgeCandidate& candidate)
    : EdgeEdgeCandidate(candidate)
{
}

EdgeEdgeFrictionConstraint::EdgeEdgeFrictionConstraint(
    const EdgeEdgeConstraint& constraint)
    : EdgeEdgeCandidate(constraint.edge0_id, constraint.edge1_id)
{
    this->weight = constraint.weight;
    this->weight_gradient = constraint.weight_gradient;
}

// ============================================================================

double EdgeEdgeFrictionConstraint::compute_distance(const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    // The distance type is known because mollified PP and PE were skipped.
    return edge_edge_distance(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()), EdgeEdgeDistanceType::EA_EB);
}

VectorMax12d EdgeEdgeFrictionConstraint::compute_distance_gradient(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    VectorMax12d grad_d;
    // The distance type is known because mollified PP and PE were skipped.
    edge_edge_distance_gradient(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()), grad_d, EdgeEdgeDistanceType::EA_EB);
    return grad_d;
}

// ============================================================================

MatrixMax<double, 3, 2>
EdgeEdgeFrictionConstraint::compute_tangent_basis(const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return edge_edge_tangent_basis(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()));
}

MatrixMax<double, 36, 2>
EdgeEdgeFrictionConstraint::compute_tangent_basis_jacobian(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return edge_edge_tangent_basis_jacobian(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()));
}

// ============================================================================

VectorMax2d
EdgeEdgeFrictionConstraint::compute_closest_point(const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return edge_edge_closest_point(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()));
}

MatrixMax<double, 2, 12>
EdgeEdgeFrictionConstraint::compute_closest_point_jacobian(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return edge_edge_closest_point_jacobian(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()));
}

// ============================================================================

MatrixMax<double, 3, 12> EdgeEdgeFrictionConstraint::relative_velocity_matrix(
    const VectorMax2d& closest_point) const
{
    assert(closest_point.size() == 2);
    return edge_edge_relative_velocity_matrix(dim(), closest_point);
}

MatrixMax<double, 6, 12>
EdgeEdgeFrictionConstraint::relative_velocity_matrix_jacobian(
    const VectorMax2d& closest_point) const
{
    assert(closest_point.size() == 2);
    return edge_edge_relative_velocity_matrix_jacobian(dim(), closest_point);
}

} // namespace ipc
