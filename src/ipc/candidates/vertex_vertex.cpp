#include "vertex_vertex.hpp"

#include <ipc/distance/point_point.hpp>

#include <iostream>

namespace ipc {

VertexVertexCandidate::VertexVertexCandidate(long _vertex0_id, long _vertex1_id)
    : vertex0_id(_vertex0_id)
    , vertex1_id(_vertex1_id)
{
}

double
VertexVertexCandidate::compute_distance(const VectorMax12d& positions) const
{
    assert(positions.size() == 4 || positions.size() == 6);
    const int dim = this->dim(positions.size());
    return point_point_distance(positions.head(dim), positions.tail(dim));
}

VectorMax12d VertexVertexCandidate::compute_distance_gradient(
    const VectorMax12d& positions) const
{
    assert(positions.size() == 4 || positions.size() == 6);
    const int dim = this->dim(positions.size());
    return point_point_distance_gradient(
        positions.head(dim), positions.tail(dim));
}

MatrixMax12d VertexVertexCandidate::compute_distance_hessian(
    const VectorMax12d& positions) const
{
    assert(positions.size() == 4 || positions.size() == 6);
    const int dim = this->dim(positions.size());
    return point_point_distance_hessian(
        positions.head(dim), positions.tail(dim));
}

bool VertexVertexCandidate::ccd(
    const VectorMax12d& vertices_t0,
    const VectorMax12d& vertices_t1,
    double& toi,
    const double min_distance,
    const double tmax,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling) const
{
    assert(vertices_t0.size() == 4 || vertices_t0.size() == 6);
    assert(vertices_t0.size() == vertices_t1.size());
    const int dim = vertices_t0.size() / 2;
    return point_point_ccd(
        vertices_t0.head(dim), vertices_t0.tail(dim), vertices_t1.head(dim),
        vertices_t1.tail(dim), toi, min_distance, tmax, tolerance,
        max_iterations, conservative_rescaling);
}

bool VertexVertexCandidate::operator==(const VertexVertexCandidate& other) const
{
    // (i, j) == (i, j) || (i, j) == (j, i)
    return (this->vertex0_id == other.vertex0_id
            && this->vertex1_id == other.vertex1_id)
        || (this->vertex0_id == other.vertex1_id
            && this->vertex1_id == other.vertex0_id);
}

bool VertexVertexCandidate::operator!=(const VertexVertexCandidate& other) const
{
    return !(*this == other);
}

bool VertexVertexCandidate::operator<(const VertexVertexCandidate& other) const
{
    long this_min = std::min(this->vertex0_id, this->vertex1_id);
    long other_min = std::min(other.vertex0_id, other.vertex1_id);
    if (this_min == other_min) {
        return std::max(this->vertex0_id, this->vertex1_id)
            < std::max(other.vertex0_id, other.vertex1_id);
    }
    return this_min < other_min;
}

} // namespace ipc
