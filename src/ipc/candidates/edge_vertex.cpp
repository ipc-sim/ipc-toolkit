#include "edge_vertex.hpp"

#include <ipc/distance/point_edge.hpp>
#include <ipc/tangent/closest_point.hpp>

#include <iostream>

namespace ipc {

EdgeVertexCandidate::EdgeVertexCandidate(index_t _edge_id, index_t _vertex_id)
    : edge_id(_edge_id)
    , vertex_id(_vertex_id)
{
}

double EdgeVertexCandidate::compute_distance(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    assert(positions.size() == 6 || positions.size() == 9);
    const int dim = this->dim(positions.size());
    return point_edge_distance(
        positions.head(dim), positions.segment(dim, dim), positions.tail(dim),
        known_dtype());
}

VectorMax12d EdgeVertexCandidate::compute_distance_gradient(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    assert(positions.size() == 6 || positions.size() == 9);
    const int dim = this->dim(positions.size());
    return point_edge_distance_gradient(
        positions.head(dim), positions.segment(dim, dim), positions.tail(dim),
        known_dtype());
}

MatrixMax12d EdgeVertexCandidate::compute_distance_hessian(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    assert(positions.size() == 6 || positions.size() == 9);
    const int dim = this->dim(positions.size());
    return point_edge_distance_hessian(
        positions.head(dim), positions.segment(dim, dim), positions.tail(dim),
        known_dtype());
}

VectorMax4d EdgeVertexCandidate::compute_coefficients(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    assert(positions.size() == 6 || positions.size() == 9);
    const int dim = this->dim(positions.size());
    Eigen::ConstRef<VectorMax3d> p = positions.head(dim);
    Eigen::ConstRef<VectorMax3d> t0 = positions.segment(dim, dim);
    Eigen::ConstRef<VectorMax3d> t1 = positions.tail(dim);

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

VectorMax3d EdgeVertexCandidate::compute_unnormalized_normal(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    const int dim = this->dim(positions.size());

    if (dim == 2) {
        // In 2D, the normal is simply the perpendicular vector to the edge
        const Eigen::Vector2d e = positions.tail<2>() - positions.segment<2>(2);
        return Eigen::Vector2d(-e.y(), e.x());
    }

    // Use triple product expansion of the cross product -e × (e × d)
    // (https://en.wikipedia.org/wiki/Cross_product#Triple_product_expansion)
    // NOTE: This would work in 2D as well, but we handle that case above.
    assert(dim == 3);
    const Eigen::Vector3d e = positions.tail<3>() - positions.segment<3>(3);
    const Eigen::Vector3d d = positions.head<3>() - positions.segment<3>(3);
    return d * e.dot(e) - e * e.dot(d);
}

MatrixMax<double, 3, 12>
EdgeVertexCandidate::compute_unnormalized_normal_jacobian(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    const int dim = this->dim(positions.size());
    if (dim == 2) {
        // In 2D, the normal is simply the perpendicular vector to the edge
        MatrixMax<double, 3, 12> dn(2, 6);
        dn.leftCols<2>().setZero();
        dn.middleCols<2>(2) << 0, 1, -1, 0;
        dn.rightCols<2>() << 0, -1, 1, 0;
        return dn;
    }

    assert(dim == 3);
    const Eigen::Vector3d e = positions.tail<3>() - positions.segment<3>(3);
    const Eigen::Vector3d d = positions.head<3>() - positions.segment<3>(3);

    const auto I = Eigen::Matrix3d::Identity();

    MatrixMax<double, 3, 12> dn(3, 9);
    // ∂n/∂x0
    dn.leftCols<3>() = e.dot(e) * I - e * e.transpose();
    // ∂n/∂x2
    dn.rightCols<3>() =
        -e.dot(d) * I - e * d.transpose() + (2 * d) * e.transpose();
    // ∂n/∂x1
    dn.middleCols<3>(3) = -dn.leftCols<3>() - dn.rightCols<3>();
    return dn;
}

bool EdgeVertexCandidate::ccd(
    Eigen::ConstRef<VectorMax12d> vertices_t0,
    Eigen::ConstRef<VectorMax12d> vertices_t1,
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
