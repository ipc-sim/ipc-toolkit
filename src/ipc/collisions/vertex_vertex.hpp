#pragma once

#include <ipc/candidates/vertex_vertex.hpp>
#include <ipc/collisions/collision_constraint.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class VertexVertexConstraint : public VertexVertexCandidate,
                               public CollisionConstraint {
public:
    using VertexVertexCandidate::VertexVertexCandidate;

    VertexVertexConstraint(const VertexVertexCandidate& candidate)
        : VertexVertexCandidate(candidate)
    {
    }

    double compute_distance(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return VertexVertexCandidate::compute_distance(vertices, edges, faces);
    }

    VectorMax12d compute_distance_gradient(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return VertexVertexCandidate::compute_distance_gradient(
            vertices, edges, faces);
    }

    MatrixMax12d compute_distance_hessian(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return VertexVertexCandidate::compute_distance_hessian(
            vertices, edges, faces);
    }

    template <typename H>
    friend H AbslHashValue(H h, const VertexVertexConstraint& vv)
    {
        return AbslHashValue(
            std::move(h), static_cast<const VertexVertexCandidate&>(vv));
    }
};

} // namespace ipc
