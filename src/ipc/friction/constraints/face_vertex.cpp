#include "face_vertex.hpp"

#include <ipc/friction/closest_point.hpp>
#include <ipc/friction/relative_velocity.hpp>
#include <ipc/friction/tangent_basis.hpp>
#include <ipc/distance/point_triangle.hpp>

namespace ipc {

FaceVertexFrictionConstraint::FaceVertexFrictionConstraint(
    long face_id, long vertex_id)
    : FaceVertexCandidate(face_id, vertex_id)
{
}

FaceVertexFrictionConstraint::FaceVertexFrictionConstraint(
    const FaceVertexCandidate& candidate)
    : FaceVertexCandidate(candidate)
{
}

FaceVertexFrictionConstraint::FaceVertexFrictionConstraint(
    const FaceVertexConstraint& constraint)
    : FaceVertexCandidate(constraint.face_id, constraint.vertex_id)
{
    this->weight = constraint.weight;
    this->weight_gradient = constraint.weight_gradient;
}

// ============================================================================

double
FaceVertexFrictionConstraint::compute_distance(const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_triangle_distance(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()), PointTriangleDistanceType::P_T);
}

VectorMax12d FaceVertexFrictionConstraint::compute_distance_gradient(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    VectorMax12d grad_d;
    point_triangle_distance_gradient(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()), grad_d, PointTriangleDistanceType::P_T);
    return grad_d;
}

// ============================================================================

MatrixMax<double, 3, 2>
FaceVertexFrictionConstraint::compute_tangent_basis(const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_triangle_tangent_basis(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()));
}

MatrixMax<double, 36, 2>
FaceVertexFrictionConstraint::compute_tangent_basis_jacobian(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_triangle_tangent_basis_jacobian(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()));
}

// ============================================================================

VectorMax2d
FaceVertexFrictionConstraint::compute_closest_point(const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_triangle_closest_point(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()));
}

MatrixMax<double, 2, 12>
FaceVertexFrictionConstraint::compute_closest_point_jacobian(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_triangle_closest_point_jacobian(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()));
}

// ============================================================================

MatrixMax<double, 3, 12> FaceVertexFrictionConstraint::relative_velocity_matrix(
    const VectorMax2d& closest_point) const
{
    assert(closest_point.size() == 2);
    return point_triangle_relative_velocity_matrix(dim(), closest_point);
}

MatrixMax<double, 6, 12>
FaceVertexFrictionConstraint::relative_velocity_matrix_jacobian(
    const VectorMax2d& closest_point) const
{
    assert(closest_point.size() == 2);
    return point_triangle_relative_velocity_matrix_jacobian(
        dim(), closest_point);
}

} // namespace ipc
