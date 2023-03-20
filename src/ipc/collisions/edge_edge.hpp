#pragma once

#include <ipc/candidates/edge_edge.hpp>
#include <ipc/collisions/collision_constraint.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class EdgeEdgeConstraint : public EdgeEdgeCandidate,
                           public CollisionConstraint {
public:
    EdgeEdgeConstraint(long edge0_id, long edge1_id, double eps_x);
    EdgeEdgeConstraint(const EdgeEdgeCandidate& candidate, double eps_x);

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
