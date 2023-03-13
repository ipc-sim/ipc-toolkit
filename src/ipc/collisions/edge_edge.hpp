#pragma once

#include <ipc/candidates/edge_edge.hpp>
#include <ipc/collisions/collision_constraint.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

struct EdgeEdgeConstraint : EdgeEdgeCandidate, CollisionConstraint {
    EdgeEdgeConstraint(long edge0_id, long edge1_id, double eps_x);
    EdgeEdgeConstraint(const EdgeEdgeCandidate& candidate, double eps_x);

    int num_vertices() const override
    {
        return EdgeEdgeCandidate::num_vertices();
    }

    std::array<long, 4> vertex_ids(
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return EdgeEdgeCandidate::vertex_ids(edges, faces);
    }

    double compute_distance(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return EdgeEdgeCandidate::compute_distance(vertices, edges, faces);
    }

    VectorMax12d compute_distance_gradient(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return EdgeEdgeCandidate::compute_distance_gradient(
            vertices, edges, faces);
    }

    MatrixMax12d compute_distance_hessian(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return EdgeEdgeCandidate::compute_distance_hessian(
            vertices, edges, faces);
    }

    double compute_potential(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat) const override;

    VectorMax12d compute_potential_gradient(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat) const override;

    MatrixMax12d compute_potential_hessian(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat,
        const bool project_hessian_to_psd) const override;

    template <typename H>
    friend H AbslHashValue(H h, const EdgeEdgeConstraint& ee)
    {
        return AbslHashValue(
            std::move(h), static_cast<const EdgeEdgeCandidate&>(ee));
    }

    double eps_x;
};

} // namespace ipc
