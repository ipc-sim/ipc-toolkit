#include "edge_edge.hpp"

#include <ipc/distance/edge_edge.hpp>

namespace ipc {

EdgeEdgeCandidate::EdgeEdgeCandidate(long _edge0_id, long _edge1_id)
    : edge0_id(_edge0_id)
    , edge1_id(_edge1_id)
{
}

double EdgeEdgeCandidate::compute_distance(const VectorMax12d& positions) const
{
    assert(positions.size() == 12);
    return edge_edge_distance(
        positions.head<3>(), positions.segment<3>(3), positions.segment<3>(6),
        positions.tail<3>(), known_dtype());
}

VectorMax12d EdgeEdgeCandidate::compute_distance_gradient(
    const VectorMax12d& positions) const
{
    assert(positions.size() == 12);
    return edge_edge_distance_gradient(
        positions.head<3>(), positions.segment<3>(3), positions.segment<3>(6),
        positions.tail<3>(), known_dtype());
}

MatrixMax12d
EdgeEdgeCandidate::compute_distance_hessian(const VectorMax12d& positions) const
{
    assert(positions.size() == 12);
    return edge_edge_distance_hessian(
        positions.head<3>(), positions.segment<3>(3), positions.segment<3>(6),
        positions.tail<3>(), known_dtype());
}

bool EdgeEdgeCandidate::ccd(
    const VectorMax12d& vertices_t0,
    const VectorMax12d& vertices_t1,
    double& toi,
    const double min_distance,
    const double tmax,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling) const
{
    assert(vertices_t0.size() == 12 && vertices_t1.size() == 12);
    return edge_edge_ccd(
        // Edge 1 at t=0
        vertices_t0.head<3>(), vertices_t0.segment<3>(3),
        // Edge 2 at t=0
        vertices_t0.segment<3>(6), vertices_t0.tail<3>(),
        // Edge 1 at t=1
        vertices_t1.head<3>(), vertices_t1.segment<3>(3),
        // Edge 2 at t=1
        vertices_t1.segment<3>(6), vertices_t1.tail<3>(), //
        toi, min_distance, tmax, tolerance, max_iterations,
        conservative_rescaling);
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
