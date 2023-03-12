#include "vertex_vertex.hpp"

#include <ipc/distance/point_point.hpp>

namespace ipc {

VertexVertexCandidate::VertexVertexCandidate(long vertex0_id, long vertex1_id)
    : vertex0_id(vertex0_id)
    , vertex1_id(vertex1_id)
{
}

double VertexVertexCandidate::compute_distance(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces) const
{
    return point_point_distance(V.row(vertex0_id), V.row(vertex1_id));
}

VectorMax6d VertexVertexCandidate::compute_distance_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces) const
{
    VectorMax6d distance_grad;
    point_point_distance_gradient(
        V.row(vertex0_id), V.row(vertex1_id), distance_grad);
    return distance_grad;
}

MatrixMax6d VertexVertexCandidate::compute_distance_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces) const
{
    MatrixMax6d distance_hess;
    point_point_distance_hessian(
        V.row(vertex0_id), V.row(vertex1_id), distance_hess);
    return distance_hess;
}

bool VertexVertexCandidate::operator==(const VertexVertexCandidate& other) const
{
    return vertex0_id == other.vertex0_id && vertex1_id == other.vertex1_id;
}

bool VertexVertexCandidate::operator!=(const VertexVertexCandidate& other) const
{
    return !(*this == other);
}

bool VertexVertexCandidate::operator<(const VertexVertexCandidate& other) const
{
    if (vertex0_id == other.vertex0_id) {
        return vertex1_id < other.vertex1_id;
    }
    return vertex0_id < other.vertex0_id;
}

} // namespace ipc
