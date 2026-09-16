#pragma once

#include <ipc/candidates/candidate_vector.hpp>
#include <ipc/candidates/stencil_adapter.hpp>
#include <ipc/candidates/stencil_mixin.hpp>
#include <ipc/distance/distance_type.hpp>

#include <Eigen/Core>

#include <array>

namespace ipc {

/// @brief A candidate for edge-vertex collision detection.
class EdgeVertexCandidate : public StencilMixin<EdgeVertexCandidate> {
    friend class StencilMixin<EdgeVertexCandidate>;

public:
    /// @brief Construct a candidate with indeterminate edge IDs.
    /// @note Keeping this trivial is what makes the type trivially copyable.
    EdgeVertexCandidate() = default;

    EdgeVertexCandidate(index_t edge_id, index_t vertex_id);

    // ------------------------------------------------------------------------
    // Stencil

    constexpr static int num_vertices() { return 3; }

    /// @brief Get the vertex IDs for the edge-vertex pair
    /// @param edges The edge connectivity matrix
    /// @param faces The face connectivity matrix
    /// @return An array of vertex IDs in the order: [vi, e0i, e1i, -1]
    std::array<index_t, 4> vertex_ids(
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const
    {
        return { { vertex_id, edges(edge_id, 0), edges(edge_id, 1), -1 } };
    }

    using StencilMixin<EdgeVertexCandidate>::compute_coefficients;
    using StencilMixin<EdgeVertexCandidate>::compute_distance;
    using StencilMixin<EdgeVertexCandidate>::compute_distance_gradient;
    using StencilMixin<EdgeVertexCandidate>::compute_distance_hessian;

    double compute_distance(
        Eigen::ConstRef<VectorMax12d> positions,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO) const;

    VectorMax12d compute_distance_gradient(
        Eigen::ConstRef<VectorMax12d> positions,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO) const;

    MatrixMax12d compute_distance_hessian(
        Eigen::ConstRef<VectorMax12d> positions,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO) const;

    VectorMax4d compute_coefficients(
        Eigen::ConstRef<VectorMax12d> positions,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO) const;

    // ------------------------------------------------------------------------

    bool
    ccd(Eigen::ConstRef<VectorMax12d> vertices_t0,
        Eigen::ConstRef<VectorMax12d> vertices_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0,
        const NarrowPhaseCCD& narrow_phase_ccd =
            DEFAULT_NARROW_PHASE_CCD) const;

    // ------------------------------------------------------------------------

    constexpr static PointEdgeDistanceType known_dtype()
    {
        return PointEdgeDistanceType::AUTO;
    }

    bool operator==(const EdgeVertexCandidate& other) const;
    bool operator!=(const EdgeVertexCandidate& other) const;
    /// @brief Compare EdgeVertexCandidates for sorting.
    bool operator<(const EdgeVertexCandidate& other) const;

    template <typename H>
    friend H AbslHashValue(H h, const EdgeVertexCandidate& ev)
    {
        return H::combine(std::move(h), ev.edge_id, ev.vertex_id);
    }

    /// @brief ID of the edge
    index_t edge_id;
    /// @brief ID of the vertex
    index_t vertex_id;

protected:
    VectorMax3d
    compute_unnormalized_normal(Eigen::ConstRef<VectorMax12d> positions) const;

    MatrixMax<double, 3, 12> compute_unnormalized_normal_jacobian(
        Eigen::ConstRef<VectorMax12d> positions) const;
};

/// @brief EdgeVertexCandidate with the runtime-polymorphic stencil interface.
///
/// The base of the corresponding collision types.
using EdgeVertexStencil =
    DTypeStencilAdapter<EdgeVertexCandidate, PointEdgeDistanceType>;

} // namespace ipc
