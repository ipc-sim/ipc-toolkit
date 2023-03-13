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
    VectorMax12d distance_grad;
    point_point_distance_gradient(
        positions.head(dim), positions.tail(dim), distance_grad);
    return distance_grad;
}

MatrixMax12d VertexVertexCandidate::compute_distance_hessian(
    const VectorMax12d& positions) const
{
    assert(positions.size() == 4 || positions.size() == 6);
    const int dim = positions.size() / 2;
    MatrixMax12d distance_hess;
    point_point_distance_hessian(
        positions.head(dim), positions.tail(dim), distance_hess);
    return distance_hess;
}

bool VertexVertexCandidate::operator==(const VertexVertexCandidate& other) const
{
    return vertex0_id == other.vertex0_id && vertex1_id == other.vertex1_id;
}

bool VertexVertexCandidate::operator!=(const VertexVertexCandidate& other) const
{
    return !(*this == other);
}

bool VertexVertexCandidate::operator<(const VertexVertexCandidate& other) const
{
    if (vertex0_id == other.vertex0_id) {
        return vertex1_id < other.vertex1_id;
    }
    return vertex0_id < other.vertex0_id;
}

} // namespace ipc
