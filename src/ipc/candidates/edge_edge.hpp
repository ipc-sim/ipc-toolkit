#pragma once

#include <ipc/candidates/collision_stencil.hpp>
#include <ipc/candidates/continuous_collision_candidate.hpp>
#include <ipc/distance/distance_type.hpp>

#include <Eigen/Core>

#include <array>

namespace ipc {

class EdgeEdgeCandidate : public ContinuousCollisionCandidate {
public:
    EdgeEdgeCandidate(long edge0_id, long edge1_id);

    int num_vertices() const override { return 4; };

    std::array<long, 4> vertex_ids(
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return { { edges(edge0_id, 0), edges(edge0_id, 1), //
                   edges(edge1_id, 0), edges(edge1_id, 1) } };
    }

    // ------------------------------------------------------------------------

    void print_ccd_query(
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override;

    // ------------------------------------------------------------------------

    bool operator==(const EdgeEdgeCandidate& other) const;
    bool operator!=(const EdgeEdgeCandidate& other) const;
    /// @brief Compare EdgeEdgeCandidates for sorting.
    bool operator<(const EdgeEdgeCandidate& other) const;

    template <typename H>
    friend H AbslHashValue(H h, const EdgeEdgeCandidate& ee)
    {
        long min_ei = std::min(ee.edge0_id, ee.edge1_id);
        long max_ei = std::max(ee.edge0_id, ee.edge1_id);
        return H::combine(std::move(h), min_ei, max_ei);
    }

    // ------------------------------------------------------------------------

    long edge0_id; ///< @brief ID of the first edge.
    long edge1_id; ///< @brief ID of the second edge.

    using CollisionStencil::compute_distance;
    using CollisionStencil::compute_distance_gradient;
    using CollisionStencil::compute_distance_hessian;
    using ContinuousCollisionCandidate::ccd;

protected:
    double compute_distance(const VectorMax12d& positions) const override;

    VectorMax12d
    compute_distance_gradient(const VectorMax12d& positions) const override;

    MatrixMax12d
    compute_distance_hessian(const VectorMax12d& positions) const override;

    bool
    ccd(const VectorMax12d& vertices_t0,
        const VectorMax12d& vertices_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0,
        const double tolerance = DEFAULT_CCD_TOLERANCE,
        const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
        const double conservative_rescaling =
            DEFAULT_CCD_CONSERVATIVE_RESCALING) const override;

    virtual EdgeEdgeDistanceType known_dtype() const
    {
        return EdgeEdgeDistanceType::AUTO;
    }
};

} // namespace ipc
