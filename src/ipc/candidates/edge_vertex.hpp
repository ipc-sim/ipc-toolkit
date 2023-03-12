#pragma once

#include <ipc/candidates/continuous_collision_candidate.hpp>
#include <ipc/distance/distance_type.hpp>

#include <Eigen/Core>

#include <array>

namespace ipc {

struct EdgeVertexCandidate : ContinuousCollisionCandidate {
    EdgeVertexCandidate(long edge_id, long vertex_id);

    int num_vertices() const { return 3; };

    std::array<long, 4>
    vertex_ids(const Eigen::MatrixXi& edges, const Eigen::MatrixXi& faces) const
    {
        return { { vertex_id, edges(edge_id, 0), edges(edge_id, 1), -1 } };
    }

    std::array<VectorMax3d, 3> vertices(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const
    {
        return { { V.row(vertex_id), V.row(edges(edge_id, 0)),
                   V.row(edges(edge_id, 1)) } };
    }

    double compute_distance(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO) const;

    VectorMax9d compute_distance_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO) const;

    MatrixMax9d compute_distance_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO) const;

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

    bool operator==(const EdgeVertexCandidate& other) const;
    bool operator!=(const EdgeVertexCandidate& other) const;
    /// @brief Compare EdgeVertexCandidates for sorting.
    bool operator<(const EdgeVertexCandidate& other) const;

    template <typename H>
    friend H AbslHashValue(H h, const EdgeVertexCandidate& ev)
    {
        return H::combine(std::move(h), ev.edge_id, ev.vertex_id);
    }

    long edge_id;   ///< @brief ID of the edge
    long vertex_id; ///< @brief ID of the vertex
};

} // namespace ipc
