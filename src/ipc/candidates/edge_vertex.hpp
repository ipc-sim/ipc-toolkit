#pragma once

#include <ipc/candidates/collision_stencil.hpp>
#include <ipc/distance/distance_type.hpp>

#include <Eigen/Core>

#include <array>

namespace ipc {

class EdgeVertexCandidate : virtual public CollisionStencil<4> {
public:
    EdgeVertexCandidate(index_t edge_id, index_t vertex_id);

    // ------------------------------------------------------------------------
    // CollisionStencil

    int num_vertices() const override { return 3; };

    std::array<index_t, 4> vertex_ids(
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const override
    {
        return { { vertex_id, edges(edge_id, 0), edges(edge_id, 1), -1 } };
    }

    using CollisionStencil<4>::compute_coefficients;
    using CollisionStencil<4>::compute_distance;
    using CollisionStencil<4>::compute_distance_gradient;
    using CollisionStencil<4>::compute_distance_hessian;

    double
    compute_distance(Eigen::ConstRef<VectorMax12d> positions) const override;

    VectorMax12d compute_distance_gradient(
        Eigen::ConstRef<VectorMax12d> positions) const override;

    MatrixMax12d compute_distance_hessian(
        Eigen::ConstRef<VectorMax12d> positions) const override;

    VectorMax4d compute_coefficients(
        Eigen::ConstRef<VectorMax12d> positions) const override;

    // ------------------------------------------------------------------------

    bool
    ccd(Eigen::ConstRef<VectorMax12d> vertices_t0,
        Eigen::ConstRef<VectorMax12d> vertices_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0,
        const NarrowPhaseCCD& narrow_phase_ccd =
            DEFAULT_NARROW_PHASE_CCD) const override;

    // ------------------------------------------------------------------------

    virtual PointEdgeDistanceType known_dtype() const
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
};

} // namespace ipc
