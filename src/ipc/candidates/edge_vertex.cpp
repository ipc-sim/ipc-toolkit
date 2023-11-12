#include "edge_vertex.hpp"

#include <ipc/distance/point_edge.hpp>
#include <ipc/utils/save_obj.hpp>

#include <iostream>

namespace ipc {

EdgeVertexCandidate::EdgeVertexCandidate(long _edge_id, long _vertex_id)
    : edge_id(_edge_id)
    , vertex_id(_vertex_id)
{
}

double
EdgeVertexCandidate::compute_distance(const VectorMax12d& positions) const
{
    assert(positions.size() == 6 || positions.size() == 9);
    const int dim = positions.size() / 3;
    return point_edge_distance(
        positions.head(dim), positions.segment(dim, dim), positions.tail(dim),
        known_dtype());
}

VectorMax12d EdgeVertexCandidate::compute_distance_gradient(
    const VectorMax12d& positions) const
{
    assert(positions.size() == 6 || positions.size() == 9);
    const int dim = positions.size() / 3;
    return point_edge_distance_gradient(
        positions.head(dim), positions.segment(dim, dim), positions.tail(dim),
        known_dtype());
}

MatrixMax12d EdgeVertexCandidate::compute_distance_hessian(
    const VectorMax12d& positions) const
{
    assert(positions.size() == 6 || positions.size() == 9);
    const int dim = positions.size() / 3;
    return point_edge_distance_hessian(
        positions.head(dim), positions.segment(dim, dim), positions.tail(dim),
        known_dtype());
}

bool EdgeVertexCandidate::ccd(
    const VectorMax12d& vertices_t0,
    const VectorMax12d& vertices_t1,
    double& toi,
    const double min_distance,
    const double tmax,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling) const
{
    assert(vertices_t0.size() == 6 || vertices_t0.size() == 9);
    assert(vertices_t0.size() == vertices_t1.size());
    const int dim = vertices_t0.size() / 3;
    return point_edge_ccd(
        // Point at t=0
        vertices_t0.head(dim),
        // Edge at t=0
        vertices_t0.segment(dim, dim), vertices_t0.tail(dim),
        // Point at t=1
        vertices_t1.head(dim),
        // Edge at t=1
        vertices_t1.segment(dim, dim), vertices_t1.tail(dim), //
        toi, min_distance, tmax, tolerance, max_iterations,
        conservative_rescaling);
}

std::ostream& EdgeVertexCandidate::write_ccd_query(
    std::ostream& out,
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces) const
{
    return out << vertices_t0.row(edges(edge_id, 0)).format(OBJ_VERTEX_FORMAT)
               << vertices_t0.row(edges(edge_id, 1)).format(OBJ_VERTEX_FORMAT)
               << vertices_t0.row(vertex_id).format(OBJ_VERTEX_FORMAT)
               << vertices_t1.row(edges(edge_id, 0)).format(OBJ_VERTEX_FORMAT)
               << vertices_t1.row(edges(edge_id, 1)).format(OBJ_VERTEX_FORMAT)
               << vertices_t1.row(vertex_id).format(OBJ_VERTEX_FORMAT);
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
