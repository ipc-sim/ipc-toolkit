#include "edge_edge.hpp"

#include <ipc/distance/edge_edge.hpp>
#include <ipc/utils/save_obj.hpp>

#include <iostream>

namespace ipc {

EdgeEdgeCandidate::EdgeEdgeCandidate(long edge0_index, long edge1_index)
    : edge0_index(edge0_index)
    , edge1_index(edge1_index)
{
}

double EdgeEdgeCandidate::compute_distance(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const EdgeEdgeDistanceType dtype) const
{
    // The distance type is unknown because of mollified PP and PE
    // constraints where also added as EE constraints.
    return edge_edge_distance(
        V.row(E(edge0_index, 0)), V.row(E(edge0_index, 1)),
        V.row(E(edge1_index, 0)), V.row(E(edge1_index, 1)), dtype);
}

VectorMax12d EdgeEdgeCandidate::compute_distance_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const EdgeEdgeDistanceType dtype) const
{
    VectorMax12d distance_grad;
    edge_edge_distance_gradient(
        V.row(E(edge0_index, 0)), V.row(E(edge0_index, 1)),
        V.row(E(edge1_index, 0)), V.row(E(edge1_index, 1)), distance_grad,
        dtype);
    return distance_grad;
}

MatrixMax12d EdgeEdgeCandidate::compute_distance_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const EdgeEdgeDistanceType dtype) const
{
    MatrixMax12d distance_hess;
    edge_edge_distance_hessian(
        V.row(E(edge0_index, 0)), V.row(E(edge0_index, 1)),
        V.row(E(edge1_index, 0)), V.row(E(edge1_index, 1)), distance_hess,
        dtype);
    return distance_hess;
}

bool EdgeEdgeCandidate::ccd(
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
    return edge_edge_ccd(
        // Edge 1 at t=0
        V0.row(E(edge0_index, 0)), V0.row(E(edge0_index, 1)),
        // Edge 2 at t=0
        V0.row(E(edge1_index, 0)), V0.row(E(edge1_index, 1)),
        // Edge 1 at t=1
        V1.row(E(edge0_index, 0)), V1.row(E(edge0_index, 1)),
        // Edge 2 at t=1
        V1.row(E(edge1_index, 0)), V1.row(E(edge1_index, 1)), //
        toi, min_distance, tmax, tolerance, max_iterations,
        conservative_rescaling);
}

void EdgeEdgeCandidate::print_ccd_query(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    std::cout << V0.row(E(edge0_index, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V0.row(E(edge0_index, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << V0.row(E(edge1_index, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V0.row(E(edge1_index, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(E(edge0_index, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(E(edge0_index, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(E(edge1_index, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(E(edge1_index, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << std::flush;
}

bool EdgeEdgeCandidate::operator==(const EdgeEdgeCandidate& other) const
{
    // (i, j) == (i, j) || (i, j) == (j, i)
    return (this->edge0_index == other.edge0_index
            && this->edge1_index == other.edge1_index)
        || (this->edge0_index == other.edge1_index
            && this->edge1_index == other.edge0_index);
}

bool EdgeEdgeCandidate::operator!=(const EdgeEdgeCandidate& other) const
{
    return !(*this == other);
}

bool EdgeEdgeCandidate::operator<(const EdgeEdgeCandidate& other) const
{
    long this_min = std::min(this->edge0_index, this->edge1_index);
    long other_min = std::min(other.edge0_index, other.edge1_index);
    if (this_min == other_min) {
        return std::max(this->edge0_index, this->edge1_index)
            < std::max(other.edge0_index, other.edge1_index);
    }
    return this_min < other_min;
}

} // namespace ipc
