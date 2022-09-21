#include "face_vertex.hpp"

#include <ipc/distance/point_triangle.hpp>

namespace ipc {

FaceVertexConstraint::FaceVertexConstraint(long face_index, long vertex_index)
    : FaceVertexCandidate(face_index, vertex_index)
{
}

FaceVertexConstraint::FaceVertexConstraint(const FaceVertexCandidate& candidate)
    : FaceVertexCandidate(candidate)
{
}

double FaceVertexConstraint::compute_distance(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    // The distance type is known because of Constraints::build()
    return point_triangle_distance(
        V.row(vertex_index), V.row(F(face_index, 0)), V.row(F(face_index, 1)),
        V.row(F(face_index, 2)), PointTriangleDistanceType::P_T);
}

VectorMax12d FaceVertexConstraint::compute_distance_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    VectorMax12d distance_grad;
    point_triangle_distance_gradient(
        V.row(vertex_index), V.row(F(face_index, 0)), V.row(F(face_index, 1)),
        V.row(F(face_index, 2)), PointTriangleDistanceType::P_T, distance_grad);
    return distance_grad;
}

MatrixMax12d FaceVertexConstraint::compute_distance_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    MatrixMax12d distance_hess;
    point_triangle_distance_hessian(
        V.row(vertex_index), V.row(F(face_index, 0)), V.row(F(face_index, 1)),
        V.row(F(face_index, 2)), PointTriangleDistanceType::P_T, distance_hess);
    return distance_hess;
}

} // namespace ipc
