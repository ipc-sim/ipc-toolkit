#include "edge_vertex.hpp"

#include <ipc/distance/point_edge.hpp>
#include <ipc/tangent/closest_point.hpp>

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
    const int dim = this->dim(positions.size());
    return point_edge_distance(
        positions.head(dim), positions.segment(dim, dim), positions.tail(dim),
        known_dtype());
}

VectorMax12d EdgeVertexCandidate::compute_distance_gradient(
    const VectorMax12d& positions) const
{
    assert(positions.size() == 6 || positions.size() == 9);
    const int dim = this->dim(positions.size());
    return point_edge_distance_gradient(
        positions.head(dim), positions.segment(dim, dim), positions.tail(dim),
        known_dtype());
}

MatrixMax12d EdgeVertexCandidate::compute_distance_hessian(
    const VectorMax12d& positions) const
{
    assert(positions.size() == 6 || positions.size() == 9);
    const int dim = this->dim(positions.size());
    return point_edge_distance_hessian(
        positions.head(dim), positions.segment(dim, dim), positions.tail(dim),
        known_dtype());
}

VectorMax4d
EdgeVertexCandidate::compute_coefficients(const VectorMax12d& positions) const
{
    assert(positions.size() == 6 || positions.size() == 9);
    const int dim = this->dim(positions.size());
    const Eigen::Ref<const VectorMax3d> p = positions.head(dim);
    const Eigen::Ref<const VectorMax3d> t0 = positions.segment(dim, dim);
    const Eigen::Ref<const VectorMax3d> t1 = positions.tail(dim);

    auto dtype = known_dtype();
    if (dtype == PointEdgeDistanceType::AUTO) {
        dtype = point_edge_distance_type(p, t0, t1);
    }

    VectorMax4d coeffs(3);
    switch (dtype) {
    case PointEdgeDistanceType::P_E0:
        coeffs << 1.0, -1.0, 0.0;
        break;
    case PointEdgeDistanceType::P_E1:
        coeffs << 1.0, 0.0, -1.0;
        break;
    case PointEdgeDistanceType::P_E: {
        const double alpha = point_edge_closest_point(p, t0, t1);
        coeffs << 1.0, -1.0 + alpha, -alpha;
    } break;
    default:
        assert(false);
        break;
    }

    return coeffs;
}

bool EdgeVertexCandidate::ccd(
    const VectorMax12d& vertices_t0,
    const VectorMax12d& vertices_t1,
    double& toi,
    const double min_distance,
    const double tmax,
    const NarrowPhaseCCD& narrow_phase_ccd) const
{
    assert(vertices_t0.size() == 6 || vertices_t0.size() == 9);
    assert(vertices_t0.size() == vertices_t1.size());
    const int dim = vertices_t0.size() / 3;
    return narrow_phase_ccd.point_edge_ccd(
        // Point at t=0
        vertices_t0.head(dim),
        // Edge at t=0
        vertices_t0.segment(dim, dim), vertices_t0.tail(dim),
        // Point at t=1
        vertices_t1.head(dim),
        // Edge at t=1
        vertices_t1.segment(dim, dim), vertices_t1.tail(dim), //
        toi, min_distance, tmax);
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
