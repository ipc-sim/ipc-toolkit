#pragma once

#include <ipc/candidates/edge_vertex.hpp>
#include <ipc/collisions/tangential/tangential_collision.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class EdgeVertexTangentialCollision : public EdgeVertexCandidate,
                                      public TangentialCollision {
public:
    using EdgeVertexCandidate::EdgeVertexCandidate;

    EdgeVertexTangentialCollision(const EdgeVertexNormalCollision& collision);

    EdgeVertexTangentialCollision(
        const EdgeVertexNormalCollision& collision,
        Eigen::ConstRef<VectorMax12d> positions,
        const NormalPotential& normal_potential,
        const double normal_stiffness);
        
    int dim() const override { return 3; }
    int ndof() const override { return 9; } // 3 vertices * 3 coordinates each 
    int num_vertices() const override { return 3; } // Edge-vertex has 3 vertices
    
    std::array<long, 4> vertex_ids(
        Eigen::ConstRef<Eigen::MatrixXi> edges, 
        Eigen::ConstRef<Eigen::MatrixXi> faces) const override
    {
        return EdgeVertexCandidate::vertex_ids(edges, faces);
    }
    
    VectorMax12d dof(
        Eigen::ConstRef<Eigen::MatrixXd> dof,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const override 
    {
        return EdgeVertexCandidate::dof(dof, edges, faces);
    }
    
    double compute_distance(Eigen::ConstRef<VectorMax12d> positions) const override
    {
        return EdgeVertexCandidate::compute_distance(positions);
    }
    
    VectorMax12d compute_distance_gradient(
        Eigen::ConstRef<VectorMax12d> positions) const override
    {
        return EdgeVertexCandidate::compute_distance_gradient(positions);
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
