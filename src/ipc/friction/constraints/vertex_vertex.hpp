#pragma once

#include <ipc/friction/constraints/friction_constraint.hpp>
#include <ipc/candidates/vertex_vertex.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class VertexVertexFrictionConstraint : public VertexVertexCandidate,
                                       public FrictionConstraint {
public:
    using VertexVertexCandidate::VertexVertexCandidate;

    VertexVertexFrictionConstraint(const VertexVertexConstraint& constraint);

    VertexVertexFrictionConstraint(
        const VertexVertexConstraint& constraint,
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

    VectorMax3d
    relative_velocity(const VectorMax12d& velocities) const override;

    using FrictionConstraint::relative_velocity_matrix;

    MatrixMax<double, 3, 12>
    relative_velocity_matrix(const VectorMax2d& closest_point) const override;

    MatrixMax<double, 6, 12> relative_velocity_matrix_jacobian(
        const VectorMax2d& closest_point) const override;
};

} // namespace ipc
