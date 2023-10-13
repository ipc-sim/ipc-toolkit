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

double
FaceVertexCandidate::compute_distance(const VectorMax12d& positions) const
{
    assert(positions.size() == 12);
    return point_triangle_distance(
        positions.head<3>(), positions.segment<3>(3), positions.segment<3>(6),
        positions.tail<3>(), known_dtype());
}

VectorMax12d FaceVertexCandidate::compute_distance_gradient(
    const VectorMax12d& positions) const
{
    assert(positions.size() == 12);
    return point_triangle_distance_gradient(
        positions.head<3>(), positions.segment<3>(3), positions.segment<3>(6),
        positions.tail<3>(), known_dtype());
}

MatrixMax12d FaceVertexCandidate::compute_distance_hessian(
    const VectorMax12d& positions) const
{
    assert(positions.size() == 12);
    return point_triangle_distance_hessian(
        positions.head<3>(), positions.segment<3>(3), positions.segment<3>(6),
        positions.tail<3>(), known_dtype());
}

bool FaceVertexCandidate::ccd(
    const VectorMax12d& vertices_t0,
    const VectorMax12d& vertices_t1,
    double& toi,
    const double min_distance,
    const double tmax,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling) const
{
    assert(vertices_t0.size() == 12 && vertices_t1.size() == 12);
    return point_triangle_ccd(
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
