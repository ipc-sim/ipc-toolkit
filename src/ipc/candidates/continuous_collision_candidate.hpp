#pragma once

#include <ipc/ccd/ccd.hpp>
#include <ipc/candidates/collision_stencil.hpp>

#include <vector>
#include <ostream>

namespace ipc {

/// Virtual class for candidates that support CCD.
class ContinuousCollisionCandidate : virtual public CollisionStencil {
public:
    virtual ~ContinuousCollisionCandidate() = default;

    /// @brief Perform narrow-phase CCD on the candidate.
    /// @param vertices_t0 Stencil vertices at the start of the time step.
    /// @param vertices_t1 Stencil vertices at the end of the time step.
    /// @param toi Computed time of impact (normalized).
    /// @param min_distance Minimum separation distance between primitives.
    /// @param tmax Maximum time (normalized) to look for collisions.
    /// @param[in] tolerance CCD tolerance used by Tight-Inclusion CCD.
    /// @param[in] max_iterations Maximum iterations used by Tight-Inclusion CCD.
    /// @param[in] conservative_rescaling Conservative rescaling value used to avoid taking steps exactly to impact.
    /// @return If the candidate had a collision over the time interval.
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

    /// @brief Write the CCD query to a stream.
    /// @param out Stream to write to.
    /// @param vertices_t0 Stencil vertices at the start of the time step.
    /// @param vertices_t1 Stencil vertices at the end of the time step.
    /// @return The stream.
    std::ostream& write_ccd_query(
        std::ostream& out,
        const VectorMax12d& vertices_t0,
        const VectorMax12d& vertices_t1) const;
};

} // namespace ipc
