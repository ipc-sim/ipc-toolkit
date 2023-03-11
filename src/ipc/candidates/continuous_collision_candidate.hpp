#pragma once

#include <ipc/ccd/ccd.hpp>

#include <Eigen/Core>

#include <vector>

namespace ipc {

/// Virtual class for candidates that support CCD.
struct ContinuousCollisionCandidate {
    virtual ~ContinuousCollisionCandidate() { }

    /// Perform narrow-phase CCD on the candidate.
    /// @param[in] V0 Mesh vertex positions at the start of the time step.
    /// @param[in] V1 Mesh vertex positions at the end of the time step.
    /// @param[in] E Mesh edges as rows of indicies into V.
    /// @param[in] F Mesh triangular faces as rows of indicies into V.
    /// @param[out] toi Computed time of impact (normalized).
    /// @param[in] min_distance Minimum separation distance between primitives.
    /// @param[in] tmax Maximum time (normalized) to look for collisions. Should be in [0, 1].
    /// @param[in] tolerance CCD tolerance used by Tight-Inclusion CCD.
    /// @param[in] max_iterations Maximum iterations used by Tight-Inclusion CCD.
    /// @param[in] conservative_rescaling Conservative rescaling value used to avoid taking steps exactly to impact.
    /// @return If the candidate had a collision over the time interval.
    virtual bool
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
            DEFAULT_CCD_CONSERVATIVE_RESCALING) const = 0;

    // Print the vertices of the CCD query for debugging.
    virtual void print_ccd_query(
        const Eigen::MatrixXd& V0,
        const Eigen::MatrixXd& V1,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const = 0;
};

} // namespace ipc
