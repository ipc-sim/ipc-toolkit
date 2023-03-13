#pragma once

#include <ipc/candidates/edge_vertex.hpp>
#include <ipc/collisions/collision_constraint.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

struct EdgeVertexConstraint : EdgeVertexCandidate, CollisionConstraint {
    using EdgeVertexCandidate::EdgeVertexCandidate;

    EdgeVertexConstraint(const EdgeVertexCandidate& candidate)
        : EdgeVertexCandidate(candidate)
    {
    }

    int num_vertices() const override { return 3; };
    std::array<long, 4> vertex_ids(
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return { { vertex_id, edges(edge_id, 0), edges(edge_id, 1), -1 } };
    }

    double compute_distance(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        // The distance type is known because of Constraints::build()
        return EdgeVertexCandidate::compute_distance(
            vertices, edges, faces, PointEdgeDistanceType::P_E);
    }

    VectorMax12d compute_distance_gradient(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        // The distance type is known because of Constraints::build()
        return EdgeVertexCandidate::compute_distance_gradient(
            vertices, edges, faces, PointEdgeDistanceType::P_E);
    }

    MatrixMax12d compute_distance_hessian(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        // The distance type is known because of Constraints::build()
        return EdgeVertexCandidate::compute_distance_hessian(
            vertices, edges, faces, PointEdgeDistanceType::P_E);
    }

    template <typename H>
    friend H AbslHashValue(H h, const EdgeVertexConstraint& ev)
    {
        return AbslHashValue(
            std::move(h), static_cast<const EdgeVertexCandidate&>(ev));
    }
};

} // namespace ipc
