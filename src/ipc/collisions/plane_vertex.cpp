#include "plane_vertex.hpp"

#include <ipc/distance/point_plane.hpp>

namespace ipc {

PlaneVertexCollision::PlaneVertexCollision(
    const VectorMax3d& _plane_origin,
    const VectorMax3d& _plane_normal,
    const long _vertex_id)
    : plane_origin(_plane_origin)
    , plane_normal(_plane_normal)
    , vertex_id(_vertex_id)
{
}

double PlaneVertexCollision::compute_distance(const VectorMax12d& point) const
{
    assert(point.size() == plane_origin.size());
    return point_plane_distance(point, plane_origin, plane_normal);
}

VectorMax12d
PlaneVertexCollision::compute_distance_gradient(const VectorMax12d& point) const
{
    assert(point.size() == plane_origin.size());
    return point_plane_distance_gradient(point, plane_origin, plane_normal);
}

MatrixMax12d
PlaneVertexCollision::compute_distance_hessian(const VectorMax12d& point) const
{
    assert(point.size() == plane_origin.size());
    return point_plane_distance_hessian(point, plane_origin, plane_normal);
}

} // namespace ipc
