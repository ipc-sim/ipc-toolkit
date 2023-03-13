#include "face_vertex.hpp"

#include <ipc/distance/point_triangle.hpp>
#include <ipc/utils/save_obj.hpp>

#include <iostream>

namespace ipc {

FaceVertexCandidate::FaceVertexCandidate(long face_id, long vertex_id)
    : face_id(face_id)
    , vertex_id(vertex_id)
{
}

double FaceVertexCandidate::compute_distance(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const PointTriangleDistanceType dtype) const
{
    return point_triangle_distance(
        vertices.row(vertex_id), vertices.row(faces(face_id, 0)),
        vertices.row(faces(face_id, 1)), vertices.row(faces(face_id, 2)),
        dtype);
}

VectorMax12d FaceVertexCandidate::compute_distance_gradient(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const PointTriangleDistanceType dtype) const
{
    VectorMax12d distance_grad;
    point_triangle_distance_gradient(
        vertices.row(vertex_id), vertices.row(faces(face_id, 0)),
        vertices.row(faces(face_id, 1)), vertices.row(faces(face_id, 2)),
        distance_grad, dtype);
    return distance_grad;
}

MatrixMax12d FaceVertexCandidate::compute_distance_hessian(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const PointTriangleDistanceType dtype) const
{
    MatrixMax12d distance_hess;
    point_triangle_distance_hessian(
        vertices.row(vertex_id), vertices.row(faces(face_id, 0)),
        vertices.row(faces(face_id, 1)), vertices.row(faces(face_id, 2)),
        distance_hess, dtype);
    return distance_hess;
}

bool FaceVertexCandidate::ccd(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    double& toi,
    const double min_distance,
    const double tmax,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling) const
{
    return point_triangle_ccd(
        // Point at t=0
        vertices_t0.row(vertex_id),
        // Triangle at t=0
        vertices_t0.row(faces(face_id, 0)), vertices_t0.row(faces(face_id, 1)),
        vertices_t0.row(faces(face_id, 2)),
        // Point at t=1
        vertices_t1.row(vertex_id),
        // Triangle at t=1
        vertices_t1.row(faces(face_id, 0)), vertices_t1.row(faces(face_id, 1)),
        vertices_t1.row(faces(face_id, 2)), //
        toi, min_distance, tmax, tolerance, max_iterations,
        conservative_rescaling);
}

void FaceVertexCandidate::print_ccd_query(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces) const
{
    std::cout << vertices_t0.row(faces(face_id, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << vertices_t0.row(faces(face_id, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << vertices_t0.row(faces(face_id, 2)).format(OBJ_VERTEX_FORMAT);
    std::cout << vertices_t0.row(vertex_id).format(OBJ_VERTEX_FORMAT);
    std::cout << vertices_t1.row(faces(face_id, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << vertices_t1.row(faces(face_id, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << vertices_t1.row(faces(face_id, 2)).format(OBJ_VERTEX_FORMAT);
    std::cout << vertices_t1.row(vertex_id).format(OBJ_VERTEX_FORMAT);
    std::cout << std::flush;
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
