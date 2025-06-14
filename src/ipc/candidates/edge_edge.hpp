#pragma once

#include <ipc/candidates/collision_stencil.hpp>
#include <ipc/distance/distance_type.hpp>

#include <Eigen/Core>

#include <array>

namespace ipc {

class EdgeEdgeCandidate : virtual public CollisionStencil<4> {
public:
    EdgeEdgeCandidate(index_t edge0_id, index_t edge1_id);

    // ------------------------------------------------------------------------
    // CollisionStencil

    int num_vertices() const override { return 4; };

    std::array<index_t, 4> vertex_ids(
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const override
    {
        return { { edges(edge0_id, 0), edges(edge0_id, 1), //
                   edges(edge1_id, 0), edges(edge1_id, 1) } };
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

    virtual EdgeEdgeDistanceType known_dtype() const
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
};

} // namespace ipc
