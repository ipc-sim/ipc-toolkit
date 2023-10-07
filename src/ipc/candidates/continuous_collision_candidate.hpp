#pragma once

#include <ipc/ccd/ccd.hpp>
#include <ipc/candidates/collision_stencil.hpp>

#include <vector>
#include <ostream>

namespace ipc {

/// Virtual class for candidates that support CCD.
class ContinuousCollisionCandidate : virtual public CollisionStencil {
public:
    virtual ~ContinuousCollisionCandidate() { }

    /// Perform narrow-phase CCD on the candidate.
    /// @param[in] vertices_t0 Mesh vertices at the start of the time step.
    /// @param[in] vertices_t1 Mesh vertices at the end of the time step.
    /// @param[in] edges Collision mesh edges as rows of indicies into vertices.
    /// @param[in] faces Collision mesh triangular faces as rows of indicies into vertices.
    /// @param[out] toi Computed time of impact (normalized).
    /// @param[in] min_distance Minimum separation distance between primitives.
    /// @param[in] tmax Maximum time (normalized) to look for collisions. Should be in [0, 1].
    /// @param[in] tolerance CCD tolerance used by Tight-Inclusion CCD.
    /// @param[in] max_iterations Maximum iterations used by Tight-Inclusion CCD.
    /// @param[in] conservative_rescaling Conservative rescaling value used to avoid taking steps exactly to impact.
    /// @return If the candidate had a collision over the time interval.
    bool
    ccd(const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0,
        const double tolerance = DEFAULT_CCD_TOLERANCE,
        const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
        const double conservative_rescaling =
            DEFAULT_CCD_CONSERVATIVE_RESCALING) const;

    /// @brief Write the CCD query to a stream.
    /// @param out Stream to write to.
    /// @param vertices_t0 Mesh vertices at the start of the time step.
    /// @param vertices_t1 Mesh vertices at the end of the time step.
    /// @param edges Collision mesh edges as rows of indicies into vertices.
    /// @param faces Collision mesh triangular faces as rows of indicies into vertices.
    /// @return The stream.
    virtual std::ostream& write_ccd_query(
        std::ostream& out,
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const = 0;

protected:
    virtual bool
    ccd(const VectorMax12d& vertices_t0,
        const VectorMax12d& vertices_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0,
        const double tolerance = DEFAULT_CCD_TOLERANCE,
        const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
        const double conservative_rescaling =
            DEFAULT_CCD_CONSERVATIVE_RESCALING) const = 0;
};

} // namespace ipc
