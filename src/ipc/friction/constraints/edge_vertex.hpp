#pragma once

#include <ipc/friction/constraints/friction_constraint.hpp>
#include <ipc/candidates/edge_vertex.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class EdgeVertexFrictionConstraint : public EdgeVertexCandidate,
                                     public FrictionConstraint {
public:
    using EdgeVertexCandidate::EdgeVertexCandidate;

    EdgeVertexFrictionConstraint(const EdgeVertexConstraint& constraint);

    EdgeVertexFrictionConstraint(
        const EdgeVertexConstraint& constraint,
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat,
        const double barrier_stiffness);

protected:
    MatrixMax<double, 3, 2>
    compute_tangent_basis(const VectorMax12d& positions) const override;

    MatrixMax<double, 36, 2> compute_tangent_basis_jacobian(
        const VectorMax12d& positions) const override;

    VectorMax2d
    compute_closest_point(const VectorMax12d& positions) const override;

    MatrixMax<double, 2, 12> compute_closest_point_jacobian(
        const VectorMax12d& positions) const override;

    VectorMax3d relative_velocity(const VectorMax12d& velocities) const override
    {
        return relative_velocity_T(velocities);
    }

    using FrictionConstraint::relative_velocity_matrix;

    MatrixMax<double, 3, 12>
    relative_velocity_matrix(const VectorMax2d& closest_point) const override;

    MatrixMax<double, 6, 12> relative_velocity_matrix_jacobian(
        const VectorMax2d& closest_point) const override;

private:
    template <typename T>
    VectorMax3<T> relative_velocity_T(const VectorMax12<T>& velocities) const
    {
        assert(velocities.size() == ndof());
        return point_edge_relative_velocity(
            velocities.head(dim()), velocities.segment(dim(), dim()),
            velocities.tail(dim()), T(closest_point[0]));
    }
};

} // namespace ipc
