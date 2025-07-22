#include "plane_vertex.hpp"

#include <ipc/ccd/point_static_plane.hpp>
#include <ipc/distance/point_plane.hpp>

namespace ipc {

PlaneVertexNormalCollision::PlaneVertexNormalCollision(
    Eigen::ConstRef<VectorMax3d> _plane_origin,
    Eigen::ConstRef<VectorMax3d> _plane_normal,
    const index_t _vertex_id)
    : plane_origin(_plane_origin)
    , plane_normal(_plane_normal)
    , vertex_id(_vertex_id)
{
}

double PlaneVertexNormalCollision::compute_distance(
    Eigen::ConstRef<VectorMax12d> point) const
{
    assert(point.size() == plane_origin.size());
    return point_plane_distance(point, plane_origin, plane_normal);
}

VectorMax12d PlaneVertexNormalCollision::compute_distance_gradient(
    Eigen::ConstRef<VectorMax12d> point) const
{
    assert(point.size() == plane_origin.size());
    return point_plane_distance_gradient(point, plane_origin, plane_normal);
}

MatrixMax12d PlaneVertexNormalCollision::compute_distance_hessian(
    Eigen::ConstRef<VectorMax12d> point) const
{
    assert(point.size() == plane_origin.size());
    return point_plane_distance_hessian(point, plane_origin, plane_normal);
}

VectorMax4d PlaneVertexNormalCollision::compute_coefficients(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    VectorMax4d coeffs(1);
    coeffs << 1.0;
    return coeffs;
}

VectorMax3d PlaneVertexNormalCollision::compute_unnormalized_normal(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    return plane_normal;
}

MatrixMax<double, 3, 12>
PlaneVertexNormalCollision::compute_unnormalized_normal_jacobian(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    return MatrixMax<double, 3, 12>::Zero(positions.size(), positions.size());
}

bool PlaneVertexNormalCollision::ccd(
    Eigen::ConstRef<VectorMax12d> vertices_t0,
    Eigen::ConstRef<VectorMax12d> vertices_t1,
    double& toi,
    const double min_distance,
    const double tmax,
    const NarrowPhaseCCD& narrow_phase_ccd) const
{
    assert(min_distance == 0 && "Not implemented");
    return point_static_plane_ccd(
        vertices_t0, vertices_t1, plane_origin, plane_normal, toi);
}

} // namespace ipc
