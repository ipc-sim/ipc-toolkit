#include "plane_vertex.hpp"

#include <ipc/ccd/point_static_plane.hpp>
#include <ipc/distance/point_plane.hpp>

namespace ipc {

PlaneVertexNormalCollision::PlaneVertexNormalCollision(
    const VectorMax3d& _plane_origin,
    const VectorMax3d& _plane_normal,
    const long _vertex_id)
    : plane_origin(_plane_origin)
    , plane_normal(_plane_normal)
    , vertex_id(_vertex_id)
{
}

double
PlaneVertexNormalCollision::compute_distance(const VectorMax12d& point) const
{
    assert(point.size() == plane_origin.size());
    return point_plane_distance(point, plane_origin, plane_normal);
}

VectorMax12d PlaneVertexNormalCollision::compute_distance_gradient(
    const VectorMax12d& point) const
{
    assert(point.size() == plane_origin.size());
    return point_plane_distance_gradient(point, plane_origin, plane_normal);
}

MatrixMax12d PlaneVertexNormalCollision::compute_distance_hessian(
    const VectorMax12d& point) const
{
    assert(point.size() == plane_origin.size());
    return point_plane_distance_hessian(point, plane_origin, plane_normal);
}

VectorMax4d PlaneVertexNormalCollision::compute_coefficients(
    const VectorMax12d& positions) const
{
    VectorMax4d coeffs(1);
    coeffs << 1.0;
    return coeffs;
}

bool PlaneVertexNormalCollision::ccd(
    const VectorMax12d& vertices_t0,
    const VectorMax12d& vertices_t1,
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
