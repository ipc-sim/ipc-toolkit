#include "plane_vertex.hpp"

#include <ipc/distance/point_triangle.hpp>
#include <ipc/tangent/closest_point.hpp>
#include <ipc/tangent/relative_velocity.hpp>
#include <ipc/tangent/tangent_basis.hpp>

namespace ipc {

PlaneVertexTangentialCollision::PlaneVertexTangentialCollision(
    const PlaneVertexNormalCollision& collision)
    : PlaneVertexCandidate(collision)
{
    this->weight = collision.weight;
    this->weight_gradient = collision.weight_gradient;
}

PlaneVertexTangentialCollision::PlaneVertexTangentialCollision(
    const PlaneVertexNormalCollision& collision,
    Eigen::ConstRef<VectorMax12d> positions,
    const double normal_force)
    : PlaneVertexTangentialCollision(collision)
{
    TangentialCollision::init(collision, positions, normal_force);
}

PlaneVertexTangentialCollision::PlaneVertexTangentialCollision(
    const PlaneVertexNormalCollision& collision,
    Eigen::ConstRef<VectorMax12d> positions,
    const NormalPotential& normal_potential)
    : PlaneVertexTangentialCollision(collision)
{
    TangentialCollision::init(collision, positions, normal_potential);
}

// ============================================================================

MatrixMax<double, 3, 2> PlaneVertexTangentialCollision::compute_tangent_basis(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    assert(positions.size() == ndof());
    Eigen::Vector3d p0 = -plane.offset() * plane.normal();
    Eigen::Vector3d p1 = p0 + plane.normal();
    return point_point_tangent_basis(
        p0.head(positions.size()), p1.head(positions.size()));
}

MatrixMax<double, 6, 12>
PlaneVertexTangentialCollision::compute_tangent_basis_jacobian(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    assert(positions.size() == ndof());
    return MatrixMax<double, 6, 12>::Zero(
        positions.size() == 2 ? 2 : 6, ndof());
}

// ============================================================================

VectorMax2d PlaneVertexTangentialCollision::compute_closest_point(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    return VectorMax2d();
}

MatrixMax<double, 2, 12>
PlaneVertexTangentialCollision::compute_closest_point_jacobian(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    return MatrixMax<double, 2, 12>();
}

// ============================================================================

VectorMax3d PlaneVertexTangentialCollision::relative_velocity(
    Eigen::ConstRef<VectorMax12d> velocities) const
{
    assert(velocities.size() <= 3);
    return velocities;
}

MatrixMax<double, 3, 12>
PlaneVertexTangentialCollision::relative_velocity_jacobian(
    Eigen::ConstRef<VectorMax2d> _closest_point) const
{
    return MatrixMax<double, 3, 12>::Identity(ndof(), ndof());
}

MatrixMax<double, 3, 24>
PlaneVertexTangentialCollision::relative_velocity_dx_dbeta(
    Eigen::ConstRef<VectorMax2d> _closest_point) const
{
    return MatrixMax<double, 3, 24>::Zero(
        dim(), _closest_point.size() * ndof());
}

} // namespace ipc