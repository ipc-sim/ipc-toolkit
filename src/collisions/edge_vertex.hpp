#pragma once

#include <ipc/broad_phase/collision_candidate.hpp>
#include <ipc/collisions/collision_constraint.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

struct EdgeVertexConstraint : EdgeVertexCandidate, CollisionConstraint {
    EdgeVertexConstraint(long edge_index, long vertex_index);
    EdgeVertexConstraint(const EdgeVertexCandidate& candidate);

    int num_vertices() const override { return 3; };
    std::array<long, 4> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const override
    {
        return { { vertex_index, E(edge_index, 0), E(edge_index, 1), -1 } };
    }

    double compute_distance(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    VectorMax12d compute_distance_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    MatrixMax12d compute_distance_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;
};

} // namespace ipc

namespace std {
template <> struct hash<ipc::EdgeVertexConstraint> {
    size_t operator()(ipc::EdgeVertexConstraint const& ev) const noexcept;
};
} // namespace std
