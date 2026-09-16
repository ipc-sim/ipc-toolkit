#pragma once

#include <ipc/candidates/stencil_adapter.hpp>
#include <ipc/candidates/stencil_mixin.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>

#include <array>

namespace ipc {

/// @brief A candidate for vertex-vertex collision detection.
class VertexVertexCandidate : public StencilMixin<VertexVertexCandidate> {
    friend class StencilMixin<VertexVertexCandidate>;

public:
    /// @brief Construct a candidate with indeterminate edge IDs.
    /// @note Keeping this trivial is what makes the type trivially copyable.
    VertexVertexCandidate() = default;

    VertexVertexCandidate(index_t vertex0_id, index_t vertex1_id);

    // ------------------------------------------------------------------------
    // Stencil

    constexpr static int num_vertices() { return 2; }

    /// @brief Get the indices of the vertices
    /// @param edges edge matrix of mesh
    /// @param faces face matrix of mesh
    /// @return List of vertex indices
    std::array<index_t, 4> vertex_ids(
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const
    {
        return { { vertex0_id, vertex1_id, -1, -1 } };
    }

    using StencilMixin<VertexVertexCandidate>::compute_coefficients;
    using StencilMixin<VertexVertexCandidate>::compute_distance;
    using StencilMixin<VertexVertexCandidate>::compute_distance_gradient;
    using StencilMixin<VertexVertexCandidate>::compute_distance_hessian;

    double compute_distance(Eigen::ConstRef<VectorMax12d> positions) const;

    VectorMax12d
    compute_distance_gradient(Eigen::ConstRef<VectorMax12d> positions) const;

    MatrixMax12d
    compute_distance_hessian(Eigen::ConstRef<VectorMax12d> positions) const;

    VectorMax4d
    compute_coefficients(Eigen::ConstRef<VectorMax12d> positions) const;

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

    bool operator==(const VertexVertexCandidate& other) const;
    bool operator!=(const VertexVertexCandidate& other) const;
    /// @brief Compare EdgeVertexCandidates for sorting.
    bool operator<(const VertexVertexCandidate& other) const;

    template <typename H>
    friend H AbslHashValue(H h, const VertexVertexCandidate& vv)
    {
        index_t min_vi = std::min(vv.vertex0_id, vv.vertex1_id);
        index_t max_vi = std::max(vv.vertex0_id, vv.vertex1_id);
        return H::combine(std::move(h), min_vi, max_vi);
    }

    /// @brief ID of the first vertex
    index_t vertex0_id;
    /// @brief ID of the second vertex
    index_t vertex1_id;

protected:
    VectorMax3d
    compute_unnormalized_normal(Eigen::ConstRef<VectorMax12d> positions) const;

    MatrixMax<double, 3, 12> compute_unnormalized_normal_jacobian(
        Eigen::ConstRef<VectorMax12d> positions) const;
};

/// @brief VertexVertexCandidate with the runtime-polymorphic stencil interface.
///
/// The base of the corresponding collision types.
using VertexVertexStencil = StencilAdapter<VertexVertexCandidate>;

} // namespace ipc
