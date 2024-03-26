#pragma once

#include <ipc/friction/collisions/friction_collision.hpp>
#include <ipc/candidates/edge_edge.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class EdgeEdgeFrictionCollision : public EdgeEdgeCandidate,
                                  public FrictionCollision {
public:
    using EdgeEdgeCandidate::EdgeEdgeCandidate;

    EdgeEdgeFrictionCollision(const EdgeEdgeCollision& collision);

    EdgeEdgeFrictionCollision(
        const EdgeEdgeCollision& collision,
        const VectorMax12d& positions,
        const BarrierPotential& barrier_potential,
        const double barrier_stiffness);

protected:
    EdgeEdgeDistanceType known_dtype() const override
    {
        // The distance type is known because mollified PP and PE were skipped.
        return EdgeEdgeDistanceType::EA_EB;
    }

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

    using FrictionCollision::relative_velocity_matrix;

    MatrixMax<double, 3, 12>
    relative_velocity_matrix(const VectorMax2d& closest_point) const override;

    MatrixMax<double, 6, 12> relative_velocity_matrix_jacobian(
        const VectorMax2d& closest_point) const override;
};

} // namespace ipc
