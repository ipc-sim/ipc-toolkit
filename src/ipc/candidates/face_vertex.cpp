#include "face_vertex.hpp"

#include <ipc/distance/point_triangle.hpp>
#include <ipc/tangent/closest_point.hpp>

#include <iostream>

namespace ipc {

FaceVertexCandidate::FaceVertexCandidate(index_t _face_id, index_t _vertex_id)
    : face_id(_face_id)
    , vertex_id(_vertex_id)
{
}

double FaceVertexCandidate::compute_distance(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    assert(positions.size() == 12);
    return point_triangle_distance(
        positions.head<3>(), positions.segment<3>(3), positions.segment<3>(6),
        positions.tail<3>(), known_dtype());
}

VectorMax12d FaceVertexCandidate::compute_distance_gradient(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    assert(positions.size() == 12);
    return point_triangle_distance_gradient(
        positions.head<3>(), positions.segment<3>(3), positions.segment<3>(6),
        positions.tail<3>(), known_dtype());
}

MatrixMax12d FaceVertexCandidate::compute_distance_hessian(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    assert(positions.size() == 12);
    return point_triangle_distance_hessian(
        positions.head<3>(), positions.segment<3>(3), positions.segment<3>(6),
        positions.tail<3>(), known_dtype());
}

VectorMax4d FaceVertexCandidate::compute_coefficients(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    assert(positions.size() == 12);
    Eigen::ConstRef<Eigen::Vector3d> p = positions.head<3>();
    Eigen::ConstRef<Eigen::Vector3d> t0 = positions.segment<3>(3);
    Eigen::ConstRef<Eigen::Vector3d> t1 = positions.segment<3>(6);
    Eigen::ConstRef<Eigen::Vector3d> t2 = positions.tail<3>();

    // Project the point inside the triangle
    auto dtype = known_dtype();
    if (dtype == PointTriangleDistanceType::AUTO) {
        dtype = point_triangle_distance_type(p, t0, t1, t2);
    }

    VectorMax4d coeffs(4);
    switch (dtype) {
    case PointTriangleDistanceType::P_T0:
        coeffs << 1.0, -1.0, 0.0, 0.0;
        break;
    case PointTriangleDistanceType::P_T1:
        coeffs << 1.0, 0.0, -1.0, 0.0;
        break;
    case PointTriangleDistanceType::P_T2:
        coeffs << 1.0, 0.0, 0.0, -1.0;
        break;
    case PointTriangleDistanceType::P_E0: {
        const double alpha = point_edge_closest_point(p, t0, t1);
        coeffs << 1.0, -1 + alpha, -alpha, 0.0;
    } break;
    case PointTriangleDistanceType::P_E1: {
        const double alpha = point_edge_closest_point(p, t1, t2);
        coeffs << 1.0, 0.0, -1 + alpha, -alpha;
    } break;
    case PointTriangleDistanceType::P_E2: {
        const double alpha = point_edge_closest_point(p, t2, t0);
        coeffs << 1.0, -alpha, 0.0, -1 + alpha;
    } break;
    case PointTriangleDistanceType::P_T: {
        const Eigen::Vector2d uv = point_triangle_closest_point(p, t0, t1, t2);
        coeffs << 1.0, -1.0 + uv[0] + uv[1], -uv[0], -uv[1];
    } break;
    default:
        assert(false);
        break;
    }

    return coeffs;
}

VectorMax3d FaceVertexCandidate::compute_unnormalized_normal(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    assert(positions.size() == 12);
    return (positions.segment<3>(6) - positions.segment<3>(3))
        .cross(positions.tail<3>() - positions.segment<3>(3));
}

MatrixMax<double, 3, 12>
FaceVertexCandidate::compute_unnormalized_normal_jacobian(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    assert(positions.size() == 12);
    MatrixMax<double, 3, 12> dn(3, 12);
    dn.leftCols<3>().setZero();
    dn.middleCols<3>(3) =
        cross_product_matrix(positions.tail<3>() - positions.segment<3>(6));
    dn.middleCols<3>(6) =
        cross_product_matrix(positions.segment<3>(3) - positions.tail<3>());
    dn.rightCols<3>() =
        cross_product_matrix(positions.segment<3>(6) - positions.segment<3>(3));
    return dn;
}

bool FaceVertexCandidate::ccd(
    Eigen::ConstRef<VectorMax12d> vertices_t0,
    Eigen::ConstRef<VectorMax12d> vertices_t1,
    double& toi,
    const double min_distance,
    const double tmax,
    const NarrowPhaseCCD& narrow_phase_ccd) const
{
    assert(vertices_t0.size() == 12 && vertices_t1.size() == 12);
    return narrow_phase_ccd.point_triangle_ccd(
        // Point at t=0
        vertices_t0.head<3>(),
        // Triangle at t=0
        vertices_t0.segment<3>(3), vertices_t0.segment<3>(6),
        vertices_t0.tail<3>(),
        // Point at t=1
        vertices_t1.head<3>(),
        // Triangle at t=1
        vertices_t1.segment<3>(3), vertices_t1.segment<3>(6),
        vertices_t1.tail<3>(), //
        toi, min_distance, tmax);
}

bool FaceVertexCandidate::operator==(const FaceVertexCandidate& other) const
{
    return face_id == other.face_id && vertex_id == other.vertex_id;
}

bool FaceVertexCandidate::operator!=(const FaceVertexCandidate& other) const
{
    return !(*this == other);
}

bool FaceVertexCandidate::operator<(const FaceVertexCandidate& other) const
{
    if (face_id == other.face_id) {
        return vertex_id < other.vertex_id;
    }
    return face_id < other.face_id;
}

} // namespace ipc
