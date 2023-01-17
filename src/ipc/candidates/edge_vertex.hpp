#pragma once

#include <ipc/candidates/continuous_collision_candidate.hpp>
#include <ipc/distance/distance_type.hpp>

#include <Eigen/Core>

#include <array>

namespace ipc {

struct EdgeVertexCandidate : ContinuousCollisionCandidate {
    EdgeVertexCandidate(long edge_index, long vertex_index);

    int num_vertices() const { return 3; };

    std::array<long, 4>
    vertex_indices(const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const
    {
        return { { vertex_index, E(edge_index, 0), E(edge_index, 1), -1 } };
    }

    double compute_distance(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO) const;

    VectorMax9d compute_distance_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO) const;

    MatrixMax9d compute_distance_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO) const;

    // ------------------------------------------------------------------------

    /// Perform narrow-phase CCD on the candidate.
    /// @param[in] V0 Mesh vertex positions at the start of the time step.
    /// @param[in] V1 Mesh vertex positions at the end of the time step.
    /// @param[in] E Mesh edges as rows of indicies into V.
    /// @param[in] F Mesh triangular faces as rows of indicies into V.
    /// @param[out] toi Computed time of impact (normalized).
    /// @param[in] tmax Maximum time (normalized) to look for collisions. Should be in [0, 1].
    /// @param[in] tolerance CCD tolerance used by Tight-Inclusion CCD.
    /// @param[in] max_iterations Maximum iterations used by Tight-Inclusion CCD.
    /// @param[in] conservative_rescaling Conservative rescaling value used to avoid taking steps exactly to impact.
    /// @return If the candidate had a collision over the time interval.
    bool
    ccd(const Eigen::MatrixXd& V0,
        const Eigen::MatrixXd& V1,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
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
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    // ------------------------------------------------------------------------

    bool operator==(const EdgeVertexCandidate& other) const;
    bool operator!=(const EdgeVertexCandidate& other) const;
    /// @brief Compare EdgeVertexCandidates for sorting.
    bool operator<(const EdgeVertexCandidate& other) const;

    template <typename H>
    friend H AbslHashValue(H h, const EdgeVertexCandidate& ev)
    {
        return H::combine(std::move(h), ev.edge_index, ev.vertex_index);
    }

    long edge_index;
    long vertex_index;
};

} // namespace ipc
