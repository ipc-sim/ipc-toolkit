#pragma once

#include <ipc/candidates/vertex_vertex.hpp>
#include <ipc/collisions/tangential/tangential_collision.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class VertexVertexTangentialCollision : public VertexVertexCandidate,
                                        public TangentialCollision {
public:
    using VertexVertexCandidate::VertexVertexCandidate;

    VertexVertexTangentialCollision(
        const VertexVertexNormalCollision& collision);

    VertexVertexTangentialCollision(
        const VertexVertexNormalCollision& collision,
        Eigen::ConstRef<VectorMax12d> positions,
        const NormalPotential& normal_potential,
        const double normal_stiffness);

protected:
    MatrixMax<double, 3, 2> compute_tangent_basis(
        Eigen::ConstRef<VectorMax12d> positions) const override;

    MatrixMax<double, 36, 2> compute_tangent_basis_jacobian(
        Eigen::ConstRef<VectorMax12d> positions) const override;

    VectorMax2d compute_closest_point(
        Eigen::ConstRef<VectorMax12d> positions) const override;

    MatrixMax<double, 2, 12> compute_closest_point_jacobian(
        Eigen::ConstRef<VectorMax12d> positions) const override;

    VectorMax3d
    relative_velocity(Eigen::ConstRef<VectorMax12d> velocities) const override;

    using TangentialCollision::relative_velocity_matrix;

    MatrixMax<double, 3, 12> relative_velocity_matrix(
        Eigen::ConstRef<VectorMax2d> closest_point) const override;

    MatrixMax<double, 6, 12> relative_velocity_matrix_jacobian(
        Eigen::ConstRef<VectorMax2d> closest_point) const override;
};

} // namespace ipc
