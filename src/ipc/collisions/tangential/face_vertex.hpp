#pragma once

#include <ipc/candidates/face_vertex.hpp>
#include <ipc/collisions/tangential/tangential_collision.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class FaceVertexTangentialCollision : public FaceVertexCandidate,
                                      public TangentialCollision {
public:
    using FaceVertexCandidate::FaceVertexCandidate;

    FaceVertexTangentialCollision(const FaceVertexNormalCollision& collision);

    FaceVertexTangentialCollision(
        const FaceVertexNormalCollision& collision,
        Eigen::ConstRef<VectorMax12d> positions,
        const NormalPotential& normal_potential,
        const double normal_stiffness);
        
    int dim() const override { return 3; }
    int ndof() const override { return 12; } // 4 vertices * 3 coordinates each
    int num_vertices() const override { return 4; } // Face-vertex has 4 vertices
    
    std::array<long, 4> vertex_ids(
        Eigen::ConstRef<Eigen::MatrixXi> edges, 
        Eigen::ConstRef<Eigen::MatrixXi> faces) const override
    {
        return FaceVertexCandidate::vertex_ids(edges, faces);
    }
    
    VectorMax12d dof(
        Eigen::ConstRef<Eigen::MatrixXd> dof,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const override 
    {
        return FaceVertexCandidate::dof(dof, edges, faces);
    }
    
    double compute_distance(Eigen::ConstRef<VectorMax12d> positions) const override
    {
        return FaceVertexCandidate::compute_distance(positions);
    }
    
    VectorMax12d compute_distance_gradient(
        Eigen::ConstRef<VectorMax12d> positions) const override
    {
        return FaceVertexCandidate::compute_distance_gradient(positions);
    }

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
