#pragma once

#include <ipc/candidates/edge_edge.hpp>
#include <ipc/collisions/tangential/tangential_collision.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class EdgeEdgeTangentialCollision : public EdgeEdgeCandidate,
                                    public TangentialCollision {
public:
    using EdgeEdgeCandidate::EdgeEdgeCandidate;

    EdgeEdgeTangentialCollision(const EdgeEdgeNormalCollision& collision);

    EdgeEdgeTangentialCollision(
        const EdgeEdgeNormalCollision& collision,
        Eigen::ConstRef<VectorMax12d> positions,
        const NormalPotential& normal_potential,
        const double normal_stiffness);
        
    int dim() const override { return 3; } // Use fixed dimension 3 since it's a 3D problem
    int ndof() const override { return 12; } // 4 vertices * 3 coordinates each
    int num_vertices() const override { return 4; } // Edge-edge has 4 vertices
    
    std::array<long, 4> vertex_ids(
        Eigen::ConstRef<Eigen::MatrixXi> edges, 
        Eigen::ConstRef<Eigen::MatrixXi> faces) const override
    {
        return EdgeEdgeCandidate::vertex_ids(edges, faces);
    }
    
    VectorMax12d dof(
        Eigen::ConstRef<Eigen::MatrixXd> dof,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const override 
    {
        return EdgeEdgeCandidate::dof(dof, edges, faces);
    }
    
    double compute_distance(Eigen::ConstRef<VectorMax12d> positions) const override
    {
        return EdgeEdgeCandidate::compute_distance(positions);
    }
    
    VectorMax12d compute_distance_gradient(
        Eigen::ConstRef<VectorMax12d> positions) const override
    {
        return EdgeEdgeCandidate::compute_distance_gradient(positions);
    }

protected:
    EdgeEdgeDistanceType known_dtype() const override
    {
        // The distance type is known because mollified PP and PE were skipped.
        return EdgeEdgeDistanceType::EA_EB;
    }

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
