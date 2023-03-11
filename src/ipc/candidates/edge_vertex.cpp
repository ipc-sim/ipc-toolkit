#include "edge_vertex.hpp"

#include <ipc/distance/point_edge.hpp>
#include <ipc/utils/save_obj.hpp>

#include <iostream>

namespace ipc {

EdgeVertexCandidate::EdgeVertexCandidate(long edge_index, long vertex_index)
    : edge_index(edge_index)
    , vertex_index(vertex_index)
{
}

double EdgeVertexCandidate::compute_distance(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const PointEdgeDistanceType dtype) const
{
    return point_edge_distance(
        V.row(vertex_index), V.row(E(edge_index, 0)), V.row(E(edge_index, 1)),
        dtype);
}

VectorMax9d EdgeVertexCandidate::compute_distance_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const PointEdgeDistanceType dtype) const
{
    VectorMax9d distance_grad;
    point_edge_distance_gradient(
        V.row(vertex_index), V.row(E(edge_index, 0)), V.row(E(edge_index, 1)),
        distance_grad, dtype);
    return distance_grad;
}

MatrixMax9d EdgeVertexCandidate::compute_distance_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const PointEdgeDistanceType dtype) const
{
    MatrixMax9d distance_hess;
    point_edge_distance_hessian(
        V.row(vertex_index), V.row(E(edge_index, 0)), V.row(E(edge_index, 1)),
        distance_hess, dtype);
    return distance_hess;
}

bool EdgeVertexCandidate::ccd(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double& toi,
    const double min_distance,
    const double tmax,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling) const
{
    return point_edge_ccd(
        // Point at t=0
        V0.row(vertex_index),
        // Edge at t=0
        V0.row(E(edge_index, 0)), V0.row(E(edge_index, 1)),
        // Point at t=1
        V1.row(vertex_index),
        // Edge at t=1
        V1.row(E(edge_index, 0)), V1.row(E(edge_index, 1)), //
        toi, min_distance, tmax, tolerance, max_iterations,
        conservative_rescaling);
}

void EdgeVertexCandidate::print_ccd_query(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    std::cout << V0.row(E(edge_index, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V0.row(E(edge_index, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << V0.row(vertex_index).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(E(edge_index, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(E(edge_index, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(vertex_index).format(OBJ_VERTEX_FORMAT);
    std::cout << std::flush;
}

bool EdgeVertexCandidate::operator==(const EdgeVertexCandidate& other) const
{
    return edge_index == other.edge_index && vertex_index == other.vertex_index;
}

bool EdgeVertexCandidate::operator!=(const EdgeVertexCandidate& other) const
{
    return !(*this == other);
}

bool EdgeVertexCandidate::operator<(const EdgeVertexCandidate& other) const
{
    if (edge_index == other.edge_index) {
        return vertex_index < other.vertex_index;
    }
    return edge_index < other.edge_index;
}

} // namespace ipc
