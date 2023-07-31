#include "vertex_vertex.hpp"

#include <ipc/distance/point_point.hpp>

namespace ipc {

VertexVertexCandidate::VertexVertexCandidate(long vertex0_id, long vertex1_id)
    : vertex0_id(vertex0_id)
    , vertex1_id(vertex1_id)
{
}

double
VertexVertexCandidate::compute_distance(const VectorMax12d& positions) const
{
    assert(positions.size() == 4 || positions.size() == 6);
    const int dim = positions.size() / 2;
    return point_point_distance(positions.head(dim), positions.tail(dim));
}

VectorMax12d VertexVertexCandidate::compute_distance_gradient(
    const VectorMax12d& positions) const
{
    assert(positions.size() == 4 || positions.size() == 6);
    const int dim = positions.size() / 2;
    return point_point_distance_gradient(
        positions.head(dim), positions.tail(dim));
}

MatrixMax12d VertexVertexCandidate::compute_distance_hessian(
    const VectorMax12d& positions) const
{
    assert(positions.size() == 4 || positions.size() == 6);
    const int dim = positions.size() / 2;
    return point_point_distance_hessian(
        positions.head(dim), positions.tail(dim));
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
