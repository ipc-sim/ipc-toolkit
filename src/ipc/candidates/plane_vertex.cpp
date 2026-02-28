#include "plane_vertex.hpp"

#include <ipc/ccd/point_static_plane.hpp>
#include <ipc/distance/point_plane.hpp>

namespace ipc {

PlaneVertexCandidate::PlaneVertexCandidate(
    const Eigen::Hyperplane<double, 3>& _plane, const index_t _vertex_id)
    : plane(_plane)
    , vertex_id(_vertex_id)
{
}

double PlaneVertexCandidate::compute_distance(
    Eigen::ConstRef<VectorMax12d> point) const
{
    assert(point.size() == 3);
    return point_plane_distance(
        point, -plane.offset() * plane.normal(), plane.normal());
}

VectorMax12d PlaneVertexCandidate::compute_distance_gradient(
    Eigen::ConstRef<VectorMax12d> point) const
{
    assert(point.size() == 3);
    return point_plane_distance_gradient(
        point, -plane.offset() * plane.normal(), plane.normal());
}

MatrixMax12d PlaneVertexCandidate::compute_distance_hessian(
    Eigen::ConstRef<VectorMax12d> point) const
{
    assert(point.size() == 3);
    return point_plane_distance_hessian(
        point, -plane.offset() * plane.normal(), plane.normal());
}

VectorMax4d PlaneVertexCandidate::compute_coefficients(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    VectorMax4d coeffs(1);
    coeffs << 1.0;
    return coeffs;
}

VectorMax3d PlaneVertexCandidate::compute_unnormalized_normal(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    return plane.normal();
}

MatrixMax<double, 3, 12>
PlaneVertexCandidate::compute_unnormalized_normal_jacobian(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    return MatrixMax<double, 3, 12>::Zero(positions.size(), positions.size());
}

bool PlaneVertexCandidate::ccd(
    Eigen::ConstRef<VectorMax12d> vertices_t0,
    Eigen::ConstRef<VectorMax12d> vertices_t1,
    double& toi,
    const double min_distance,
    const double tmax,
    const NarrowPhaseCCD& narrow_phase_ccd) const
{
    assert(min_distance == 0 && "Not implemented");
    assert(vertices_t0.size() == 3 && vertices_t1.size() == 3);
    return point_static_plane_ccd(
        vertices_t0, vertices_t1, -plane.offset() * plane.normal(),
        plane.normal(), toi);
}

bool PlaneVertexCandidate::operator==(const PlaneVertexCandidate& other) const
{
    return this->vertex_id == other.vertex_id
        && this->plane.normal().isApprox(other.plane.normal())
        && this->plane.offset() == other.plane.offset();
}

bool PlaneVertexCandidate::operator!=(const PlaneVertexCandidate& other) const
{
    return !(*this == other);
}

bool PlaneVertexCandidate::operator<(const PlaneVertexCandidate& other) const
{
    if (this->vertex_id == other.vertex_id) {
        return this->plane.offset() < other.plane.offset();
    }
    return this->vertex_id < other.vertex_id;
}

} // namespace ipc
