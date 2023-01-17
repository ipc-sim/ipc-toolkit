#include "face_vertex.hpp"

#include <ipc/distance/point_triangle.hpp>
#include <ipc/utils/save_obj.hpp>

#include <iostream>

namespace ipc {

FaceVertexCandidate::FaceVertexCandidate(long face_index, long vertex_index)
    : face_index(face_index)
    , vertex_index(vertex_index)
{
}

double FaceVertexCandidate::compute_distance(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const PointTriangleDistanceType dtype) const
{
    return point_triangle_distance(
        V.row(vertex_index), V.row(F(face_index, 0)), V.row(F(face_index, 1)),
        V.row(F(face_index, 2)), dtype);
}

VectorMax12d FaceVertexCandidate::compute_distance_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const PointTriangleDistanceType dtype) const
{
    VectorMax12d distance_grad;
    point_triangle_distance_gradient(
        V.row(vertex_index), V.row(F(face_index, 0)), V.row(F(face_index, 1)),
        V.row(F(face_index, 2)), distance_grad, dtype);
    return distance_grad;
}

MatrixMax12d FaceVertexCandidate::compute_distance_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const PointTriangleDistanceType dtype) const
{
    MatrixMax12d distance_hess;
    point_triangle_distance_hessian(
        V.row(vertex_index), V.row(F(face_index, 0)), V.row(F(face_index, 1)),
        V.row(F(face_index, 2)), distance_hess, dtype);
    return distance_hess;
}

bool FaceVertexCandidate::ccd(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double& toi,
    const double min_distance,
    const double tmax,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling) const
{
    return point_triangle_ccd(
        // Point at t=0
        V0.row(vertex_index),
        // Triangle at t=0
        V0.row(F(face_index, 0)), V0.row(F(face_index, 1)),
        V0.row(F(face_index, 2)),
        // Point at t=1
        V1.row(vertex_index),
        // Triangle at t=1
        V1.row(F(face_index, 0)), V1.row(F(face_index, 1)),
        V1.row(F(face_index, 2)), //
        toi, min_distance, tmax, tolerance, max_iterations,
        conservative_rescaling);
}

void FaceVertexCandidate::print_ccd_query(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    std::cout << V0.row(F(face_index, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V0.row(F(face_index, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << V0.row(F(face_index, 2)).format(OBJ_VERTEX_FORMAT);
    std::cout << V0.row(vertex_index).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(F(face_index, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(F(face_index, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(F(face_index, 2)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(vertex_index).format(OBJ_VERTEX_FORMAT);
    std::cout << std::flush;
}

bool FaceVertexCandidate::operator==(const FaceVertexCandidate& other) const
{
    return face_index == other.face_index && vertex_index == other.vertex_index;
}

bool FaceVertexCandidate::operator!=(const FaceVertexCandidate& other) const
{
    return !(*this == other);
}

bool FaceVertexCandidate::operator<(const FaceVertexCandidate& other) const
{
    if (face_index == other.face_index) {
        return vertex_index < other.vertex_index;
    }
    return face_index < other.face_index;
}

} // namespace ipc
