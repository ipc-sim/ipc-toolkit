#pragma once

#include <ipc/broad_phase/collision_candidate.hpp>
#include <ipc/collisions/collision_constraint.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

struct VertexVertexConstraint : VertexVertexCandidate, CollisionConstraint {
    VertexVertexConstraint(long vertex0_index, long vertex1_index);
    VertexVertexConstraint(const VertexVertexCandidate& candidate);

    int num_vertices() const override { return 2; };
    std::array<long, 4> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const override
    {
        return { { vertex0_index, vertex1_index, -1, -1 } };
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
template <> struct hash<ipc::VertexVertexConstraint> {
    size_t operator()(ipc::VertexVertexConstraint const& vv) const noexcept;
};
} // namespace std
