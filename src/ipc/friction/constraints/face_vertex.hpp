#pragma once

#include <ipc/friction/constraints/friction_constraint.hpp>
#include <ipc/candidates/face_vertex.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

struct FaceVertexFrictionConstraint : FaceVertexCandidate, FrictionConstraint {
    FaceVertexFrictionConstraint(long face_id, long vertex_id);
    FaceVertexFrictionConstraint(const FaceVertexCandidate& constraint);
    FaceVertexFrictionConstraint(const FaceVertexConstraint& constraint);
    FaceVertexFrictionConstraint(
        const FaceVertexConstraint& constraint,
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat,
        const double barrier_stiffness)
        : FaceVertexFrictionConstraint(constraint)
    {
        FrictionConstraint::init(
            vertices, edges, faces, dhat, barrier_stiffness,
            constraint.minimum_distance);
    }

    int num_vertices() const override { return 4; }
    std::array<long, 4> vertex_ids(
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return { { vertex_id, //
                   faces(face_id, 0), faces(face_id, 1), faces(face_id, 2) } };
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
        return point_triangle_relative_velocity(
            u.head(dim()), u.segment(dim(), dim()), u.segment(2 * dim(), dim()),
            u.tail(dim()), closest_point.cast<T>());
    }
};

} // namespace ipc
