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

    std::array<long, 4> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const override
    {
        return FaceVertexCandidate::vertex_indices(E, F);
    }

    double compute_distance(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override
    {
        // The distance type is known because of Constraints::build()
        return FaceVertexCandidate::compute_distance(
            V, E, F, PointTriangleDistanceType::P_T);
    }

    VectorMax12d compute_distance_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override
    {
        // The distance type is known because of Constraints::build()
        return FaceVertexCandidate::compute_distance_gradient(
            V, E, F, PointTriangleDistanceType::P_T);
    }

    MatrixMax12d compute_distance_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override
    {
        // The distance type is known because of Constraints::build()
        return FaceVertexCandidate::compute_distance_hessian(
            V, E, F, PointTriangleDistanceType::P_T);
    }

    template <typename H>
    friend H AbslHashValue(H h, const FaceVertexConstraint& fv)
    {
        return AbslHashValue(
            std::move(h), static_cast<const FaceVertexCandidate&>(fv));
    }
};

} // namespace ipc
