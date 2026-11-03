#pragma once

#include <ipc/candidates/candidate_vector.hpp>
#include <ipc/candidates/stencil_adapter.hpp>
#include <ipc/candidates/stencil_mixin.hpp>
#include <ipc/collisions/normal/normal_collision.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class PlaneVertexCandidate : public StencilMixin<PlaneVertexCandidate> {
    friend class StencilMixin<PlaneVertexCandidate>;

public:
    /// @brief Construct a candidate with indeterminate edge IDs.
    /// @note Keeping this trivial is what makes the type trivially copyable.
    PlaneVertexCandidate() = default;

    PlaneVertexCandidate(
        const Eigen::Hyperplane<double, 3>& plane, const index_t vertex_id);

    constexpr static int num_vertices() { return 1; }

    std::array<index_t, 4> vertex_ids(
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const
    {
        return { { vertex_id, -1, -1, -1 } };
    }

    using StencilMixin<PlaneVertexCandidate>::compute_coefficients;
    using StencilMixin<PlaneVertexCandidate>::compute_distance;
    using StencilMixin<PlaneVertexCandidate>::compute_distance_gradient;
    using StencilMixin<PlaneVertexCandidate>::compute_distance_hessian;

    /// @brief Compute the distance between the point and plane.
    /// @param point Point's position.
    /// @return Distance of the stencil.
    double compute_distance(Eigen::ConstRef<VectorMax12d> point) const;

    /// @brief Compute the gradient of the distance w.r.t. the point's positions.
    /// @param point Point's position.
    /// @return Distance gradient w.r.t. the point's positions.
    VectorMax12d
    compute_distance_gradient(Eigen::ConstRef<VectorMax12d> point) const;

    /// @brief Compute the distance Hessian of the stencil w.r.t. the stencil's vertex positions.
    /// @param point Point's position.
    /// @return Distance Hessian w.r.t. the point's positions.
    MatrixMax12d
    compute_distance_hessian(Eigen::ConstRef<VectorMax12d> point) const;

    /// @brief Compute the coefficients of the stencil.
    /// @param positions Vertex positions.
    /// @return Coefficients of the stencil.
    VectorMax4d
    compute_coefficients(Eigen::ConstRef<VectorMax12d> positions) const;

    /// @brief Perform narrow-phase CCD on the candidate.
    /// @param[in] vertices_t0 Stencil vertices at the start of the time step.
    /// @param[in] vertices_t1 Stencil vertices at the end of the time step.
    /// @param[out] toi Computed time of impact (normalized).
    /// @param[in] min_distance Minimum separation distance between primitives.
    /// @param[in] tmax Maximum time (normalized) to look for collisions.
    /// @param[in] narrow_phase_ccd The narrow phase CCD algorithm to use.
    /// @return If the candidate had a collision over the time interval.
    bool
    ccd(Eigen::ConstRef<VectorMax12d> vertices_t0,
        Eigen::ConstRef<VectorMax12d> vertices_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0,
        const NarrowPhaseCCD& narrow_phase_ccd =
            DEFAULT_NARROW_PHASE_CCD) const;

    /// @brief The plane of the candidate.
    Eigen::Hyperplane<double, 3> plane;

    /// @brief The vertex's id.
    index_t vertex_id;

    /// @brief Compare two PlaneVertexCandidates for equality.
    bool operator==(const PlaneVertexCandidate& other) const;

    /// @brief Compare two PlaneVertexCandidates for inequality.
    bool operator!=(const PlaneVertexCandidate& other) const;

    /// @brief Compare two PlaneVertexCandidates for less than.
    bool operator<(const PlaneVertexCandidate& other) const;

protected:
    /// @brief Compute the normal vector of the stencil.
    /// @param positions Vertex positions.
    /// @return Normal vector of the stencil.
    VectorMax3d
    compute_unnormalized_normal(Eigen::ConstRef<VectorMax12d> positions) const;

    /// @brief Compute the Jacobian of the normal vector of the stencil.
    /// @param positions Vertex positions.
    /// @return Jacobian of the normal vector of the stencil.
    MatrixMax<double, 3, 12> compute_unnormalized_normal_jacobian(
        Eigen::ConstRef<VectorMax12d> positions) const;
};

/// @brief PlaneVertexCandidate with the runtime-polymorphic stencil interface.
///
/// The base of the corresponding collision types.
using PlaneVertexStencil = StencilAdapter<PlaneVertexCandidate>;

} // namespace ipc