#pragma once

#include <ipc/friction/constraints/friction_constraint.hpp>
#include <ipc/candidates/edge_edge.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class EdgeEdgeFrictionConstraint : public EdgeEdgeCandidate,
                                   public FrictionConstraint {
public:
    EdgeEdgeFrictionConstraint(long edge0_id, long edge1_id);
    EdgeEdgeFrictionConstraint(const EdgeEdgeCandidate& constraint);
    EdgeEdgeFrictionConstraint(const EdgeEdgeConstraint& constraint);
    EdgeEdgeFrictionConstraint(
        const EdgeEdgeConstraint& constraint,
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat,
        const double barrier_stiffness)
        : EdgeEdgeFrictionConstraint(constraint)
    {
        FrictionConstraint::init(
            vertices, edges, faces, dhat, barrier_stiffness,
            constraint.minimum_distance);
    }

    template <typename T>
    T compute_potential(
        const MatrixX<T>& velocities,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double epsv_times_h) const
    {
        return compute_potential_common(
            relative_velocity_T(select_dof(velocities, edges, faces)),
            epsv_times_h);
    }

protected:
    virtual double compute_distance(const VectorMax12d& x) const override;
    virtual VectorMax12d
    compute_distance_gradient(const VectorMax12d& x) const override;

    MatrixMax<double, 3, 2>
    compute_tangent_basis(const VectorMax12d& x) const override;

    MatrixMax<double, 36, 2>
    compute_tangent_basis_jacobian(const VectorMax12d& x) const override;

    VectorMax2d compute_closest_point(const VectorMax12d& x) const override;

    MatrixMax<double, 2, 12>
    compute_closest_point_jacobian(const VectorMax12d& x) const override;

    VectorMax3d relative_velocity(const VectorMax12d& u) const override
    {
        return relative_velocity_T(u);
    }

    using FrictionConstraint::relative_velocity_matrix;

    MatrixMax<double, 3, 12>
    relative_velocity_matrix(const VectorMax2d& closest_point) const override;

    MatrixMax<double, 6, 12> relative_velocity_matrix_jacobian(
        const VectorMax2d& closest_point) const override;

private:
    template <typename T>
    VectorMax3<T> relative_velocity_T(const VectorMax12<T>& u) const
    {
        assert(u.size() == ndof());
        return edge_edge_relative_velocity(
            u.head(dim()), u.segment(dim(), dim()), u.segment(2 * dim(), dim()),
            u.tail(dim()), closest_point.cast<T>());
    }
};

} // namespace ipc
