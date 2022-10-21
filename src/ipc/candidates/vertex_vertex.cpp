#include "vertex_vertex.hpp"

#include <ipc/distance/point_point.hpp>

namespace ipc {

VertexVertexCandidate::VertexVertexCandidate(
    long vertex0_index, long vertex1_index)
    : vertex0_index(vertex0_index)
    , vertex1_index(vertex1_index)
{
}

double VertexVertexCandidate::compute_distance(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    return point_point_distance(V.row(vertex0_index), V.row(vertex1_index));
}

VectorMax6d VertexVertexCandidate::compute_distance_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    VectorMax6d distance_grad;
    point_point_distance_gradient(
        V.row(vertex0_index), V.row(vertex1_index), distance_grad);
    return distance_grad;
}

MatrixMax6d VertexVertexCandidate::compute_distance_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    MatrixMax6d distance_hess;
    point_point_distance_hessian(
        V.row(vertex0_index), V.row(vertex1_index), distance_hess);
    return distance_hess;
}

bool VertexVertexCandidate::operator==(const VertexVertexCandidate& other) const
{
    return vertex0_index == other.vertex0_index
        && vertex1_index == other.vertex1_index;
}

bool VertexVertexCandidate::operator!=(const VertexVertexCandidate& other) const
{
    return !(*this == other);
}

bool VertexVertexCandidate::operator<(const VertexVertexCandidate& other) const
{
    if (vertex0_index == other.vertex0_index) {
        return vertex1_index < other.vertex1_index;
    }
    return vertex0_index < other.vertex0_index;
}

} // namespace ipc
