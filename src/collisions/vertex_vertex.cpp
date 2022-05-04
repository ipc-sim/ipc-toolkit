#include <ipc/collisions/vertex_vertex.hpp>

#include <ipc/distance/point_point.hpp>

namespace ipc {

VertexVertexConstraint::VertexVertexConstraint(
    long vertex0_index, long vertex1_index)
    : VertexVertexCandidate(vertex0_index, vertex1_index)
{
}

VertexVertexConstraint::VertexVertexConstraint(
    const VertexVertexCandidate& candidate)
    : VertexVertexCandidate(candidate)
{
}

double VertexVertexConstraint::compute_distance(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    return point_point_distance(V.row(vertex0_index), V.row(vertex1_index));
}

VectorMax12d VertexVertexConstraint::compute_distance_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    VectorMax6d distance_grad;
    point_point_distance_gradient(
        V.row(vertex0_index), V.row(vertex1_index), distance_grad);
    return distance_grad;
}

MatrixMax12d VertexVertexConstraint::compute_distance_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    MatrixMax6d distance_hess;
    point_point_distance_hessian(
        V.row(vertex0_index), V.row(vertex1_index), distance_hess);
    return distance_hess;
}

} // namespace ipc

namespace std {

size_t hash<ipc::VertexVertexConstraint>::operator()(
    ipc::VertexVertexConstraint const& vv) const noexcept
{
    long min_vi = std::min(vv.vertex0_index, vv.vertex1_index);
    long max_vi = std::max(vv.vertex0_index, vv.vertex1_index);
    return std::hash<long>()(min_vi) ^ std::hash<long>()(max_vi);
}

} // namespace std
