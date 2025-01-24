#include "edge_edge.hpp"

#include <ipc/distance/edge_edge.hpp>
#include <ipc/tangent/closest_point.hpp>

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

VectorMax4d
EdgeEdgeCandidate::compute_coefficients(const VectorMax12d& positions) const
{
    assert(positions.size() == 12);
    const Eigen::Ref<const Eigen::Vector3d> ea0 = positions.head<3>();
    const Eigen::Ref<const Eigen::Vector3d> ea1 = positions.segment<3>(3);
    const Eigen::Ref<const Eigen::Vector3d> eb0 = positions.segment<3>(6);
    const Eigen::Ref<const Eigen::Vector3d> eb1 = positions.tail<3>();

    // Project the point inside the triangle
    auto dtype = known_dtype();
    if (dtype == EdgeEdgeDistanceType::AUTO) {
        dtype = edge_edge_distance_type(ea0, ea1, eb0, eb1);
    }

    VectorMax4d coeffs(4);
    switch (dtype) {
    case EdgeEdgeDistanceType::EA0_EB0:
        coeffs << 1.0, 0.0, -1.0, 0.0;
        break;
    case EdgeEdgeDistanceType::EA0_EB1:
        coeffs << 1.0, 0.0, 0.0, -1.0;
        break;
    case EdgeEdgeDistanceType::EA1_EB0:
        coeffs << 0.0, 1.0, -1.0, 0.0;
        break;
    case EdgeEdgeDistanceType::EA1_EB1:
        coeffs << 0.0, 1.0, 0.0, -1.0;
        break;
    case EdgeEdgeDistanceType::EA_EB0: {
        const double alpha = point_edge_closest_point(eb0, ea0, ea1);
        coeffs << 1.0 - alpha, alpha, -1.0, 0;
    } break;
    case EdgeEdgeDistanceType::EA_EB1: {
        const double alpha = point_edge_closest_point(eb1, ea0, ea1);
        coeffs << 1.0 - alpha, alpha, 0, -1.0;
    } break;
    case EdgeEdgeDistanceType::EA0_EB: {
        const double alpha = point_edge_closest_point(ea0, eb0, eb1);
        coeffs << 1.0, 0, -1.0 + alpha, -alpha;
    } break;
    case EdgeEdgeDistanceType::EA1_EB: {
        const double alpha = point_edge_closest_point(ea1, eb0, eb1);
        coeffs << 0, 1.0, -1.0 + alpha, -alpha;
    } break;
    case EdgeEdgeDistanceType::EA_EB: {
        const Eigen::Vector2d ab = edge_edge_closest_point(ea0, ea1, eb0, eb1);
        coeffs << 1.0 - ab[0], ab[0], -1.0 + ab[1], -ab[1];
    } break;
    default:
        assert(false);
        break;
    }

    return coeffs;
}

bool EdgeEdgeCandidate::ccd(
    const VectorMax12d& vertices_t0,
    const VectorMax12d& vertices_t1,
    double& toi,
    const double min_distance,
    const double tmax,
    const NarrowPhaseCCD& narrow_phase_ccd) const
{
    assert(vertices_t0.size() == 12 && vertices_t1.size() == 12);
    return narrow_phase_ccd.edge_edge_ccd(
        // Edge 1 at t=0
        vertices_t0.head<3>(), vertices_t0.segment<3>(3),
        // Edge 2 at t=0
        vertices_t0.segment<3>(6), vertices_t0.tail<3>(),
        // Edge 1 at t=1
        vertices_t1.head<3>(), vertices_t1.segment<3>(3),
        // Edge 2 at t=1
        vertices_t1.segment<3>(6), vertices_t1.tail<3>(), //
        toi, min_distance, tmax);
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
