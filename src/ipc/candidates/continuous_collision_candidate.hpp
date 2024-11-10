#pragma once

#include <ipc/candidates/collision_stencil.hpp>
#include <ipc/ccd/default_narrow_phase_ccd.hpp>

#include <ostream>
#include <vector>

namespace ipc {

/// Virtual class for candidates that support CCD.
class ContinuousCollisionCandidate : virtual public CollisionStencil {
public:
    virtual ~ContinuousCollisionCandidate() = default;

    /// @brief Perform narrow-phase CCD on the candidate.
    /// @param[in] vertices_t0 Stencil vertices at the start of the time step.
    /// @param[in] vertices_t1 Stencil vertices at the end of the time step.
    /// @param[out] toi Computed time of impact (normalized).
    /// @param[in] min_distance Minimum separation distance between primitives.
    /// @param[in] tmax Maximum time (normalized) to look for collisions.
    /// @param[in] narrow_phase_ccd The narrow phase CCD algorithm to use.
    /// @return If the candidate had a collision over the time interval.
    virtual bool
    ccd(const VectorMax12d& vertices_t0,
        const VectorMax12d& vertices_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0,
        const NarrowPhaseCCD& narrow_phase_ccd =
            DEFAULT_NARROW_PHASE_CCD) const = 0;

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
