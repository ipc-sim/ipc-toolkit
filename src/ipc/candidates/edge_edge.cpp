#include "edge_edge.hpp"

#include <ipc/distance/edge_edge.hpp>
#include <ipc/utils/save_obj.hpp>

#include <iostream>

namespace ipc {

EdgeEdgeCandidate::EdgeEdgeCandidate(long edge0_id, long edge1_id)
    : edge0_id(edge0_id)
    , edge1_id(edge1_id)
{
}

double EdgeEdgeCandidate::compute_distance(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const EdgeEdgeDistanceType dtype) const
{
    // The distance type is unknown because of mollified PP and PE
    // constraints where also added as EE constraints.
    return edge_edge_distance(
        V.row(edges(edge0_id, 0)), V.row(edges(edge0_id, 1)),
        V.row(edges(edge1_id, 0)), V.row(edges(edge1_id, 1)), dtype);
}

VectorMax12d EdgeEdgeCandidate::compute_distance_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const EdgeEdgeDistanceType dtype) const
{
    VectorMax12d distance_grad;
    edge_edge_distance_gradient(
        V.row(edges(edge0_id, 0)), V.row(edges(edge0_id, 1)),
        V.row(edges(edge1_id, 0)), V.row(edges(edge1_id, 1)), distance_grad,
        dtype);
    return distance_grad;
}

MatrixMax12d EdgeEdgeCandidate::compute_distance_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const EdgeEdgeDistanceType dtype) const
{
    MatrixMax12d distance_hess;
    edge_edge_distance_hessian(
        V.row(edges(edge0_id, 0)), V.row(edges(edge0_id, 1)),
        V.row(edges(edge1_id, 0)), V.row(edges(edge1_id, 1)), distance_hess,
        dtype);
    return distance_hess;
}

bool EdgeEdgeCandidate::ccd(
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
    return edge_edge_ccd(
        // Edge 1 at t=0
        V0.row(edges(edge0_id, 0)), V0.row(edges(edge0_id, 1)),
        // Edge 2 at t=0
        V0.row(edges(edge1_id, 0)), V0.row(edges(edge1_id, 1)),
        // Edge 1 at t=1
        V1.row(edges(edge0_id, 0)), V1.row(edges(edge0_id, 1)),
        // Edge 2 at t=1
        V1.row(edges(edge1_id, 0)), V1.row(edges(edge1_id, 1)), //
        toi, min_distance, tmax, tolerance, max_iterations,
        conservative_rescaling);
}

void EdgeEdgeCandidate::print_ccd_query(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces) const
{
    std::cout << V0.row(edges(edge0_id, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V0.row(edges(edge0_id, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << V0.row(edges(edge1_id, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V0.row(edges(edge1_id, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(edges(edge0_id, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(edges(edge0_id, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(edges(edge1_id, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(edges(edge1_id, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << std::flush;
}

bool EdgeEdgeCandidate::operator==(const EdgeEdgeCandidate& other) const
{
    // (i, j) == (i, j) || (i, j) == (j, i)
    return (this->edge0_id == other.edge0_id
            && this->edge1_id == other.edge1_id)
        || (this->edge0_id == other.edge1_id
            && this->edge1_id == other.edge0_id);
}

bool EdgeEdgeCandidate::operator!=(const EdgeEdgeCandidate& other) const
{
    return !(*this == other);
}

bool EdgeEdgeCandidate::operator<(const EdgeEdgeCandidate& other) const
{
    long this_min = std::min(this->edge0_id, this->edge1_id);
    long other_min = std::min(other.edge0_id, other.edge1_id);
    if (this_min == other_min) {
        return std::max(this->edge0_id, this->edge1_id)
            < std::max(other.edge0_id, other.edge1_id);
    }
    return this_min < other_min;
}

} // namespace ipc
