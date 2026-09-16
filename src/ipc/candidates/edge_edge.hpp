#pragma once

#include <ipc/candidates/candidate_vector.hpp>
#include <ipc/candidates/stencil_adapter.hpp>
#include <ipc/candidates/stencil_mixin.hpp>
#include <ipc/distance/distance_type.hpp>

#include <Eigen/Core>

#include <array>

namespace ipc {

/// @brief A candidate for edge-edge collision detection.
class EdgeEdgeCandidate : public StencilMixin<EdgeEdgeCandidate> {
    friend class StencilMixin<EdgeEdgeCandidate>;

public:
    /// @brief Construct a candidate with indeterminate edge IDs.
    /// @note Keeping this trivial is what makes the type trivially copyable.
    EdgeEdgeCandidate() = default;

    EdgeEdgeCandidate(index_t edge0_id, index_t edge1_id);

    // ------------------------------------------------------------------------
    // Stencil

    constexpr static int num_vertices() { return 4; }

    /// @brief Get the vertex IDs for the edge-edge pair
    /// @param edges The edge connectivity matrix
    /// @param faces The face connectivity matrix
    /// @return An array of vertex IDs in the order: [ea0i, ea1i, eb0i, eb1i]
    std::array<index_t, 4> vertex_ids(
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const
    {
        return { { edges(edge0_id, 0), edges(edge0_id, 1), //
                   edges(edge1_id, 0), edges(edge1_id, 1) } };
    }

    using StencilMixin<EdgeEdgeCandidate>::compute_coefficients;
    using StencilMixin<EdgeEdgeCandidate>::compute_distance;
    using StencilMixin<EdgeEdgeCandidate>::compute_distance_gradient;
    using StencilMixin<EdgeEdgeCandidate>::compute_distance_hessian;

    /// @param dtype The distance type, when it is known a priori.
    /// Defaults to AUTO, which classifies the pair from the positions.
    double compute_distance(
        Eigen::ConstRef<VectorMax12d> positions,
        const EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO) const;

    VectorMax12d compute_distance_gradient(
        Eigen::ConstRef<VectorMax12d> positions,
        const EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO) const;

    MatrixMax12d compute_distance_hessian(
        Eigen::ConstRef<VectorMax12d> positions,
        const EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO) const;

    VectorMax4d compute_coefficients(
        Eigen::ConstRef<VectorMax12d> positions,
        const EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO) const;

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

    /// @brief The distance type known a priori; AUTO for a plain candidate.
    constexpr static EdgeEdgeDistanceType known_dtype()
    {
        return EdgeEdgeDistanceType::AUTO;
    }

    bool operator==(const EdgeEdgeCandidate& other) const;
    bool operator!=(const EdgeEdgeCandidate& other) const;
    /// @brief Compare EdgeEdgeCandidates for sorting.
    bool operator<(const EdgeEdgeCandidate& other) const;

    template <typename H>
    friend H AbslHashValue(H h, const EdgeEdgeCandidate& ee)
    {
        index_t min_ei = std::min(ee.edge0_id, ee.edge1_id);
        index_t max_ei = std::max(ee.edge0_id, ee.edge1_id);
        return H::combine(std::move(h), min_ei, max_ei);
    }

    /// @brief ID of the first edge.
    index_t edge0_id;
    /// @brief ID of the second edge.
    index_t edge1_id;

protected:
    VectorMax3d
    compute_unnormalized_normal(Eigen::ConstRef<VectorMax12d> positions) const;

    MatrixMax<double, 3, 12> compute_unnormalized_normal_jacobian(
        Eigen::ConstRef<VectorMax12d> positions) const;
};

/// @brief EdgeEdgeCandidate with the runtime-polymorphic stencil interface.
///
/// The base of the edge-edge collision types, which pin the distance type by
/// overriding known_dtype().
using EdgeEdgeStencil =
    DTypeStencilAdapter<EdgeEdgeCandidate, EdgeEdgeDistanceType>;

} // namespace ipc
