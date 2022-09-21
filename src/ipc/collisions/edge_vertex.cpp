#include "edge_vertex.hpp"

#include <ipc/distance/point_edge.hpp>

namespace ipc {

EdgeVertexConstraint::EdgeVertexConstraint(long edge_index, long vertex_index)
    : EdgeVertexCandidate(edge_index, vertex_index)
{
}

EdgeVertexConstraint::EdgeVertexConstraint(const EdgeVertexCandidate& candidate)
    : EdgeVertexCandidate(candidate)
{
}

double EdgeVertexConstraint::compute_distance(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    // The distance type is known because of Constraints::build()
    return point_edge_distance(
        V.row(vertex_index), V.row(E(edge_index, 0)), V.row(E(edge_index, 1)),
        PointEdgeDistanceType::P_E);
}

VectorMax12d EdgeVertexConstraint::compute_distance_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    VectorMax9d distance_grad;
    point_edge_distance_gradient(
        V.row(vertex_index), V.row(E(edge_index, 0)), V.row(E(edge_index, 1)),
        PointEdgeDistanceType::P_E, distance_grad);
    return distance_grad;
}

MatrixMax12d EdgeVertexConstraint::compute_distance_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    MatrixMax9d distance_hess;
    point_edge_distance_hessian(
        V.row(vertex_index), V.row(E(edge_index, 0)), V.row(E(edge_index, 1)),
        PointEdgeDistanceType::P_E, distance_hess);
    return distance_hess;
}

} // namespace ipc
