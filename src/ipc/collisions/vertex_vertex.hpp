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

    std::array<long, 4> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const override
    {
        return VertexVertexCandidate::vertex_indices(E, F);
    }

    double compute_distance(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override
    {
        return VertexVertexCandidate::compute_distance(V, E, F);
    }

    VectorMax12d compute_distance_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override
    {
        return VertexVertexCandidate::compute_distance_gradient(V, E, F);
    }

    MatrixMax12d compute_distance_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override
    {
        return VertexVertexCandidate::compute_distance_hessian(V, E, F);
    }

    template <typename H>
    friend H AbslHashValue(H h, const VertexVertexConstraint& vv)
    {
        return AbslHashValue(
            std::move(h), static_cast<const VertexVertexCandidate&>(vv));
    }
};

} // namespace ipc
