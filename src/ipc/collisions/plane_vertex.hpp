#pragma once

#include <ipc/collisions/collision_constraint.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class PlaneVertexConstraint : public CollisionConstraint {
public:
    PlaneVertexConstraint(
        const VectorMax3d& plane_origin,
        const VectorMax3d& plane_normal,
        const long vertex_id);

    int num_vertices() const override { return 1; };

    std::array<long, 4> vertex_ids(
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return { { vertex_id, -1, -1, -1 } };
    }

    using CollisionConstraint::compute_distance;
    using CollisionConstraint::compute_distance_gradient;
    using CollisionConstraint::compute_distance_hessian;

    VectorMax3d plane_origin;
    VectorMax3d plane_normal;
    long vertex_id;

protected:
    /// @brief Compute the distance between the point and plane.
    /// @param point Point's position.
    /// @return Distance of the stencil.
    double compute_distance(const VectorMax12d& point) const override;

    /// @brief Compute the gradient of the distance w.r.t. the point's positions.
    /// @param point Point's position.
    /// @return Distance gradient w.r.t. the point's positions.
    VectorMax12d
    compute_distance_gradient(const VectorMax12d& point) const override;

    /// @brief Compute the distance Hessian of the stencil w.r.t. the stencil's vertex positions.
    /// @param point Point's position.
    /// @return Distance Hessian w.r.t. the point's positions.
    MatrixMax12d
    compute_distance_hessian(const VectorMax12d& point) const override;
};

} // namespace ipc
