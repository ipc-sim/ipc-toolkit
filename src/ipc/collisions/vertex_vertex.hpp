#pragma once

#include <ipc/candidates/vertex_vertex.hpp>
#include <ipc/collisions/collision_constraint.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

struct VertexVertexConstraint : VertexVertexCandidate, CollisionConstraint {
    using VertexVertexCandidate::VertexVertexCandidate;

    VertexVertexConstraint(const VertexVertexCandidate& candidate)
        : VertexVertexCandidate(candidate)
    {
    }

    int num_vertices() const override
    {
        return VertexVertexCandidate::num_vertices();
    }

    std::array<long, 4> vertex_ids(
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return VertexVertexCandidate::vertex_ids(edges, faces);
    }

    double compute_distance(
        const Eigen::MatrixXd& positions,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return VertexVertexCandidate::compute_distance(positions, edges, faces);
    }

    VectorMax12d compute_distance_gradient(
        const Eigen::MatrixXd& positions,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return VertexVertexCandidate::compute_distance_gradient(
            positions, edges, faces);
    }

    MatrixMax12d compute_distance_hessian(
        const Eigen::MatrixXd& positions,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return VertexVertexCandidate::compute_distance_hessian(
            positions, edges, faces);
    }

    template <typename H>
    friend H AbslHashValue(H h, const VertexVertexConstraint& vv)
    {
        return AbslHashValue(
            std::move(h), static_cast<const VertexVertexCandidate&>(vv));
    }
};

} // namespace ipc
