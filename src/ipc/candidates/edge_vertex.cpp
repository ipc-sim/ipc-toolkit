#include "edge_vertex.hpp"

#include <ipc/distance/point_edge.hpp>
#include <ipc/utils/save_obj.hpp>

#include <iostream>

namespace ipc {

EdgeVertexCandidate::EdgeVertexCandidate(long edge_id, long vertex_id)
    : edge_id(edge_id)
    , vertex_id(vertex_id)
{
}

double EdgeVertexCandidate::compute_distance(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const PointEdgeDistanceType dtype) const
{
    return point_edge_distance(
        V.row(vertex_id), V.row(edges(edge_id, 0)), V.row(edges(edge_id, 1)),
        dtype);
}

VectorMax9d EdgeVertexCandidate::compute_distance_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const PointEdgeDistanceType dtype) const
{
    VectorMax9d distance_grad;
    point_edge_distance_gradient(
        V.row(vertex_id), V.row(edges(edge_id, 0)), V.row(edges(edge_id, 1)),
        distance_grad, dtype);
    return distance_grad;
}

MatrixMax9d EdgeVertexCandidate::compute_distance_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const PointEdgeDistanceType dtype) const
{
    MatrixMax9d distance_hess;
    point_edge_distance_hessian(
        V.row(vertex_id), V.row(edges(edge_id, 0)), V.row(edges(edge_id, 1)),
        distance_hess, dtype);
    return distance_hess;
}

bool EdgeVertexCandidate::ccd(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    double& toi,
    const double min_distance,
    const double tmax,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling) const
{
    return point_edge_ccd(
        // Point at t=0
        V0.row(vertex_id),
        // Edge at t=0
        V0.row(edges(edge_id, 0)), V0.row(edges(edge_id, 1)),
        // Point at t=1
        V1.row(vertex_id),
        // Edge at t=1
        V1.row(edges(edge_id, 0)), V1.row(edges(edge_id, 1)), //
        toi, min_distance, tmax, tolerance, max_iterations,
        conservative_rescaling);
}

void EdgeVertexCandidate::print_ccd_query(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces) const
{
    std::cout << V0.row(edges(edge_id, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V0.row(edges(edge_id, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << V0.row(vertex_id).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(edges(edge_id, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(edges(edge_id, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(vertex_id).format(OBJ_VERTEX_FORMAT);
    std::cout << std::flush;
}

bool EdgeVertexCandidate::operator==(const EdgeVertexCandidate& other) const
{
    return edge_id == other.edge_id && vertex_id == other.vertex_id;
}

bool EdgeVertexCandidate::operator!=(const EdgeVertexCandidate& other) const
{
    return !(*this == other);
}

bool EdgeVertexCandidate::operator<(const EdgeVertexCandidate& other) const
{
    if (edge_id == other.edge_id) {
        return vertex_id < other.vertex_id;
    }
    return edge_id < other.edge_id;
}

} // namespace ipc
