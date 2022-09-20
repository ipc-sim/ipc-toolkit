#include "plane_vertex.hpp"

#include <ipc/distance/point_plane.hpp>

namespace ipc {

PlaneVertexConstraint::PlaneVertexConstraint(
    const VectorMax3d& plane_origin,
    const VectorMax3d& plane_normal,
    const long vertex_index)
    : plane_origin(plane_origin)
    , plane_normal(plane_normal)
    , vertex_index(vertex_index)
{
}

double PlaneVertexConstraint::compute_distance(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    return point_plane_distance(
        V.row(vertex_index).transpose(), plane_origin, plane_normal);
}

VectorMax12d PlaneVertexConstraint::compute_distance_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    VectorMax3d distance_grad;
    point_plane_distance_gradient(
        V.row(vertex_index).transpose(), plane_origin, plane_normal,
        distance_grad);
    return distance_grad;
}

MatrixMax12d PlaneVertexConstraint::compute_distance_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    MatrixMax3d distance_hess;
    point_plane_distance_hessian(
        V.row(vertex_index).transpose(), plane_origin, plane_normal,
        distance_hess);
    return distance_hess;
}

} // namespace ipc
