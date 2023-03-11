#pragma once

#include <ipc/candidates/edge_edge.hpp>
#include <ipc/collisions/collision_constraint.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

struct EdgeEdgeConstraint : EdgeEdgeCandidate, CollisionConstraint {
    EdgeEdgeConstraint(long edge0_index, long edge1_index, double eps_x);
    EdgeEdgeConstraint(const EdgeEdgeCandidate& candidate, double eps_x);

    int num_vertices() const override
    {
        return EdgeEdgeCandidate::num_vertices();
    }

    std::array<long, 4> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const override
    {
        return EdgeEdgeCandidate::vertex_indices(E, F);
    }

    double compute_distance(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override
    {
        return EdgeEdgeCandidate::compute_distance(V, E, F);
    }

    VectorMax12d compute_distance_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override
    {
        return EdgeEdgeCandidate::compute_distance_gradient(V, E, F);
    }

    MatrixMax12d compute_distance_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override
    {
        return EdgeEdgeCandidate::compute_distance_hessian(V, E, F);
    }

    double compute_potential(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat) const override;

    VectorMax12d compute_potential_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat) const override;

    MatrixMax12d compute_potential_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
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
