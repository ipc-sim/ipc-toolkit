#pragma once

#include <ipc/friction/constraints/friction_constraint.hpp>
#include <ipc/candidates/vertex_vertex.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

struct VertexVertexFrictionConstraint : VertexVertexCandidate,
                                        FrictionConstraint {
    VertexVertexFrictionConstraint(long vertex0_id, long vertex1_id);
    VertexVertexFrictionConstraint(const VertexVertexCandidate& candidate);
    VertexVertexFrictionConstraint(const VertexVertexConstraint& constraint);
    VertexVertexFrictionConstraint(
        const VertexVertexConstraint& constraint,
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat,
        const double barrier_stiffness)
        : VertexVertexFrictionConstraint(constraint)
    {
        FrictionConstraint::init(
            V, edges, faces, dhat, barrier_stiffness,
            constraint.minimum_distance);
    }

    int num_vertices() const override
    {
        return VertexVertexCandidate::num_vertices();
    }

    std::array<long, 4> vertex_ids(
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return VertexVertexCandidate::vertex_ids(edges, faces);
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
        return point_point_relative_velocity(u.head(dim()), u.tail(dim()));
    }
};

} // namespace ipc
