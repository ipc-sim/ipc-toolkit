#pragma once

#include <ipc/candidates/candidates.hpp>
#include <ipc/ccd/nonlinear_ccd.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>

namespace ipc::rigid {

/// @brief A class for storing and managing collision candidates.
class RigidCandidates : public Candidates {
public:
    RigidCandidates() = default;

    using Candidates::build;

    /// @brief Initialize the set of discrete collision detection candidates.
    /// @param bodies The rigid bodies.
    /// @param poses The poses of the rigid bodies.
    /// @param inflation_radius Amount to inflate the bounding boxes.
    /// @param broad_phase Broad phase method to use.
    void build(
        const RigidBodies& bodies,
        const std::vector<Pose>& poses,
        const double inflation_radius = 0,
        BroadPhase* broad_phase = nullptr);

    /// @brief Initialize the set of continuous collision detection candidates.
    /// @note Assumes the trajectory is linear.
    /// @param bodies The rigid bodies.
    /// @param poses_t0 The starting poses of the rigid bodies.
    /// @param poses_t1 The ending poses of the rigid bodies.
    /// @param inflation_radius Amount to inflate the bounding boxes.
    /// @param broad_phase Broad phase method to use.
    void build(
        const RigidBodies& bodies,
        const std::vector<Pose>& poses_t0,
        const std::vector<Pose>& poses_t1,
        const double inflation_radius = 0,
        BroadPhase* broad_phase = nullptr);

    /// @brief Determine if the step is collision free from the set of candidates.
    /// @note Uses the nonlinear rigid trajectories of the bodies.
    /// @param bodies The rigid bodies.
    /// @param poses_t0 The starting poses of the rigid bodies.
    /// @param poses_t1 The ending poses of the rigid bodies.
    /// @param min_distance The minimum distance allowable between any two elements.
    /// @param nonlinear_ccd The nonlinear narrow phase CCD algorithm to use.
    /// @returns True if <b>no</b> collisions occur.
    bool is_step_collision_free(
        const RigidBodies& bodies,
        const std::vector<Pose>& poses_t0,
        const std::vector<Pose>& poses_t1,
        const double min_distance = 0.0,
        const NonlinearCCD& nonlinear_ccd = NonlinearCCD()) const;

    /// @brief Computes a maximal step size that is collision free using the set of collision candidates.
    /// @note Uses the nonlinear rigid trajectories of the bodies.
    /// @param bodies The rigid bodies.
    /// @param poses_t0 The starting poses of the rigid bodies. Assumed to be intersection free.
    /// @param poses_t1 The ending poses of the rigid bodies.
    /// @param min_distance The minimum distance allowable between any two elements.
    /// @param nonlinear_ccd The nonlinear narrow phase CCD algorithm to use.
    /// @returns A step-size \f$\in [0, 1]\f$ that is collision free. A value of 1.0 if a full step and 0.0 is no step.
    double compute_collision_free_stepsize(
        const RigidBodies& bodies,
        const std::vector<Pose>& poses_t0,
        const std::vector<Pose>& poses_t1,
        const double min_distance = 0.0,
        const NonlinearCCD& nonlinear_ccd = NonlinearCCD()) const;

    /// @brief Computes a conservative bound on the largest-feasible step size for surface primitives not in collision.
    ///
    /// Each vertex's displacement is bounded by the arc-length bound
    /// ‖Δp‖ + ‖Δθ‖‖x̄‖ of its rigid trajectory.
    ///
    /// @param bodies The rigid bodies.
    /// @param poses_t0 The starting poses of the rigid bodies.
    /// @param poses_t1 The ending poses of the rigid bodies.
    /// @param dhat Barrier activation distance.
    double compute_noncandidate_conservative_stepsize(
        const RigidBodies& bodies,
        const std::vector<Pose>& poses_t0,
        const std::vector<Pose>& poses_t1,
        const double dhat) const;

    /// @brief Computes a CFL-inspired CCD maximum step step size.
    /// @param bodies The rigid bodies.
    /// @param poses_t0 The starting poses of the rigid bodies.
    /// @param poses_t1 The ending poses of the rigid bodies.
    /// @param dhat Barrier activation distance.
    /// @param min_distance The minimum distance allowable between any two elements.
    /// @param broad_phase The broad phase algorithm to use.
    /// @param nonlinear_ccd The nonlinear narrow phase CCD algorithm to use.
    double compute_cfl_stepsize(
        const RigidBodies& bodies,
        const std::vector<Pose>& poses_t0,
        const std::vector<Pose>& poses_t1,
        const double dhat,
        const double min_distance = 0.0,
        BroadPhase* broad_phase = nullptr,
        const NonlinearCCD& nonlinear_ccd = NonlinearCCD()) const;

private:
    /// @brief Perform narrow-phase nonlinear CCD on the i-th candidate.
    /// @param[in] i Index of the candidate (same ordering as operator[]).
    /// @param[in] bodies The rigid bodies.
    /// @param[in] poses_t0 The starting poses of the rigid bodies.
    /// @param[in] poses_t1 The ending poses of the rigid bodies.
    /// @param[in] min_distance The minimum distance allowable between any two elements.
    /// @param[in] tmax Maximum time (normalized) to look for collisions.
    /// @param[in] nonlinear_ccd The nonlinear narrow phase CCD algorithm to use.
    /// @param[out] toi Computed time of impact (normalized).
    /// @return If the candidate had a collision over the time interval.
    bool candidate_ccd(
        const size_t i,
        const RigidBodies& bodies,
        const std::vector<Pose>& poses_t0,
        const std::vector<Pose>& poses_t1,
        const double min_distance,
        const double tmax,
        const NonlinearCCD& nonlinear_ccd,
        double& toi) const;
};

} // namespace ipc::rigid
