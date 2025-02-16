#pragma once

#include <ipc/collisions/normal/normal_collision.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class PlaneVertexNormalCollision : public NormalCollision {
public:
    PlaneVertexNormalCollision(
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

    using CollisionStencil::compute_coefficients;
    using CollisionStencil::compute_distance;
    using CollisionStencil::compute_distance_gradient;
    using CollisionStencil::compute_distance_hessian;

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

    /// @brief Compute the coefficients of the stencil.
    /// @param positions Vertex positions.
    /// @return Coefficients of the stencil.
    VectorMax4d
    compute_coefficients(const VectorMax12d& positions) const override;

    /// @brief Perform narrow-phase CCD on the candidate.
    /// @param[in] vertices_t0 Stencil vertices at the start of the time step.
    /// @param[in] vertices_t1 Stencil vertices at the end of the time step.
    /// @param[out] toi Computed time of impact (normalized).
    /// @param[in] min_distance Minimum separation distance between primitives.
    /// @param[in] tmax Maximum time (normalized) to look for collisions.
    /// @param[in] narrow_phase_ccd The narrow phase CCD algorithm to use.
    /// @return If the candidate had a collision over the time interval.
    bool
    ccd(const VectorMax12d& vertices_t0,
        const VectorMax12d& vertices_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0,
        const NarrowPhaseCCD& narrow_phase_ccd =
            DEFAULT_NARROW_PHASE_CCD) const override;

    /// @brief The plane's origin.
    VectorMax3d plane_origin;

    /// @brief The plane's normal.
    VectorMax3d plane_normal;

    /// @brief The vertex's id.
    long vertex_id;
};

} // namespace ipc
