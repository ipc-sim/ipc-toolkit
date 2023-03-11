#include "plane_vertex.hpp"

#include <ipc/distance/point_plane.hpp>

namespace ipc {

PlaneVertexConstraint::PlaneVertexConstraint(
    const VectorMax3d& plane_origin,
    const VectorMax3d& plane_normal,
    const long vertex_id)
    : plane_origin(plane_origin)
    , plane_normal(plane_normal)
    , vertex_id(vertex_id)
{
}

double PlaneVertexConstraint::compute_distance(
    const Eigen::MatrixXd& positions,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces) const
{
    return point_plane_distance(
        positions.row(vertex_id).transpose(), plane_origin, plane_normal);
}

VectorMax12d PlaneVertexConstraint::compute_distance_gradient(
    const Eigen::MatrixXd& positions,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces) const
{
    VectorMax3d distance_grad;
    point_plane_distance_gradient(
        positions.row(vertex_id).transpose(), plane_origin, plane_normal,
        distance_grad);
    return distance_grad;
}

MatrixMax12d PlaneVertexConstraint::compute_distance_hessian(
    const Eigen::MatrixXd& positions,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces) const
{
    MatrixMax3d distance_hess;
    point_plane_distance_hessian(
        positions.row(vertex_id).transpose(), plane_origin, plane_normal,
        distance_hess);
    return distance_hess;
}

} // namespace ipc
