#pragma once

#include <ipc/candidates/face_vertex.hpp>
#include <ipc/collisions/collision_constraint.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

struct FaceVertexConstraint : FaceVertexCandidate, CollisionConstraint {
    using FaceVertexCandidate::FaceVertexCandidate;

    FaceVertexConstraint(const FaceVertexCandidate& candidate)
        : FaceVertexCandidate(candidate)
    {
    }

    int num_vertices() const override
    {
        return FaceVertexCandidate::num_vertices();
    };

    std::array<long, 4> vertex_ids(
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return FaceVertexCandidate::vertex_ids(edges, faces);
    }

    double compute_distance(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        // The distance type is known because of Constraints::build()
        return FaceVertexCandidate::compute_distance(
            vertices, edges, faces, PointTriangleDistanceType::P_T);
    }

    VectorMax12d compute_distance_gradient(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        // The distance type is known because of Constraints::build()
        return FaceVertexCandidate::compute_distance_gradient(
            vertices, edges, faces, PointTriangleDistanceType::P_T);
    }

    MatrixMax12d compute_distance_hessian(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        // The distance type is known because of Constraints::build()
        return FaceVertexCandidate::compute_distance_hessian(
            vertices, edges, faces, PointTriangleDistanceType::P_T);
    }

    template <typename H>
    friend H AbslHashValue(H h, const FaceVertexConstraint& fv)
    {
        return AbslHashValue(
            std::move(h), static_cast<const FaceVertexCandidate&>(fv));
    }
};

} // namespace ipc
