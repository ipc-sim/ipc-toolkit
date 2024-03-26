#pragma once

#include <ipc/collisions/collision.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class PlaneVertexCollision : public Collision {
public:
    PlaneVertexCollision(
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

    /// @brief The plane's origin.
    VectorMax3d plane_origin;

    /// @brief The plane's normal.
    VectorMax3d plane_normal;

    /// @brief The vertex's id.
    long vertex_id;
};

} // namespace ipc
