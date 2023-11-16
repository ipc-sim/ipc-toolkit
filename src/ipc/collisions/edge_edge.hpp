#pragma once

#include <ipc/candidates/edge_edge.hpp>
#include <ipc/collisions/collision_constraint.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class EdgeEdgeConstraint : public EdgeEdgeCandidate,
                           public CollisionConstraint {
public:
    EdgeEdgeConstraint(
        const long edge0_id,
        const long edge1_id,
        const double eps_x,
        const EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO);

    EdgeEdgeConstraint(
        const EdgeEdgeCandidate& candidate,
        const double eps_x,
        const EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO);

    EdgeEdgeConstraint(
        const long edge0_id,
        const long edge1_id,
        const double eps_x,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient,
        const EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO);

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

    // ------------------------------------------------------------------------

    bool operator==(const EdgeEdgeConstraint& other) const;
    bool operator!=(const EdgeEdgeConstraint& other) const;
    bool operator<(const EdgeEdgeConstraint& other) const;

    template <typename H>
    friend H AbslHashValue(H h, const EdgeEdgeConstraint& ee)
    {
        return H::combine(
            std::move(h), static_cast<const EdgeEdgeCandidate&>(ee), ee.dtype);
    }

    // ------------------------------------------------------------------------

    /// @brief Mollifier activation threshold.
    /// @see edge_edge_mollifier
    double eps_x;

    /// @brief Cached distance type.
    /// Some EE constraints are mollified EV or VV constraints.
    EdgeEdgeDistanceType dtype;

protected:
    virtual EdgeEdgeDistanceType known_dtype() const override { return dtype; }

    MatrixMax12d compute_shape_derivative_second_term(
        const Eigen::MatrixXd& rest_positions,
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat) const override;
};

} // namespace ipc
