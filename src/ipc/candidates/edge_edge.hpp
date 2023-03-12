#pragma once

#include <ipc/candidates/continuous_collision_candidate.hpp>
#include <ipc/distance/distance_type.hpp>

#include <Eigen/Core>

#include <array>

namespace ipc {

struct EdgeEdgeCandidate : ContinuousCollisionCandidate {
    EdgeEdgeCandidate(long edge0_id, long edge1_id);

    int num_vertices() const { return 4; };

    std::array<long, 4>
    vertex_ids(const Eigen::MatrixXi& edges, const Eigen::MatrixXi& faces) const
    {
        return { { edges(edge0_id, 0), edges(edge0_id, 1), //
                   edges(edge1_id, 0), edges(edge1_id, 1) } };
    }

    std::array<Eigen::Vector3d, 4> vertices(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const
    {
        assert(V.cols() == 3);
        return { {
            V.row(edges(edge0_id, 0)),
            V.row(edges(edge0_id, 1)),
            V.row(edges(edge1_id, 0)),
            V.row(edges(edge1_id, 1)),
        } };
    }

    double compute_distance(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO) const;

    VectorMax12d compute_distance_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO) const;

    MatrixMax12d compute_distance_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO) const;

    // ------------------------------------------------------------------------

    /// Perform narrow-phase CCD on the candidate.
    /// @param[in] V0 Mesh vertex V at the start of the time step.
    /// @param[in] V1 Mesh vertex V at the end of the time step.
    /// @param[in] edges Mesh edges as rows of indicies into V.
    /// @param[in] faces Mesh triangular faces as rows of indicies into V.
    /// @param[out] toi Computed time of impact (normalized).
    /// @param[in] tmax Maximum time (normalized) to look for collisions. Should be in [0, 1].
    /// @param[in] tolerance CCD tolerance used by Tight-Inclusion CCD.
    /// @param[in] max_iterations Maximum iterations used by Tight-Inclusion CCD.
    /// @param[in] conservative_rescaling Conservative rescaling value used to avoid taking steps exactly to impact.
    /// @return If the candidate had a collision over the time interval.
    bool
    ccd(const Eigen::MatrixXd& V0,
        const Eigen::MatrixXd& V1,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0,
        const double tolerance = DEFAULT_CCD_TOLERANCE,
        const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
        const double conservative_rescaling =
            DEFAULT_CCD_CONSERVATIVE_RESCALING) const override;

    void print_ccd_query(
        const Eigen::MatrixXd& V0,
        const Eigen::MatrixXd& V1,
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
};

} // namespace ipc
