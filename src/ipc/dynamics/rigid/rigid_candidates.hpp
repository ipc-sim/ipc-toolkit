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
    /// @note Assumes the trajectory is linear.
    /// @param bodies The rigid bodies.
    /// @param poses_t0 The starting poses of the rigid bodies.
    /// @param poses_t1 The ending poses of the rigid bodies.
    /// @param min_distance The minimum distance allowable between any two elements.
    /// @param narrow_phase_ccd The narrow phase CCD algorithm to use.
    /// @returns True if <b>any</b> collisions occur.
    bool is_step_collision_free(
        const CollisionMesh& mesh,
        const std::vector<Pose>& poses_t0,
        const std::vector<Pose>& poses_t1,
        const double min_distance = 0.0,
        const NonlinearCCD& nonlinear_ccd = NonlinearCCD()) const;

    /// @brief Computes a maximal step size that is collision free using the set of collision candidates.
    /// @note Assumes the trajectory is linear.
    /// @param bodies The rigid bodies.
    /// @param poses_t0 The starting poses of the rigid bodies.
    /// @param poses_t1 The ending poses of the rigid bodies.
    /// @param min_distance The minimum distance allowable between any two elements.
    /// @param narrow_phase_ccd The narrow phase CCD algorithm to use.
    /// @returns A step-size \f$\in [0, 1]\f$ that is collision free. A value of 1.0 if a full step and 0.0 is no step.
    double compute_collision_free_stepsize(
        const CollisionMesh& mesh,
        const std::vector<Pose>& poses_t0,
        const std::vector<Pose>& poses_t1,
        const double min_distance = 0.0,
        const NonlinearCCD& nonlinear_ccd = NonlinearCCD()) const;

    /// @brief Computes a conservative bound on the largest-feasible step size for surface primitives not in collision.
    /// @param bodies The rigid bodies.
    /// @param displacements Surface vertex displacements (rowwise).
    /// @param dhat Barrier activation distance.
    double compute_noncandidate_conservative_stepsize(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> displacements,
        const double dhat) const;

    /// @brief Computes a CFL-inspired CCD maximum step step size.
    /// @param bodies The rigid bodies.
    /// @param poses_t0 The starting poses of the rigid bodies.
    /// @param poses_t1 The ending poses of the rigid bodies.
    /// @param dhat Barrier activation distance.
    /// @param min_distance The minimum distance allowable between any two elements.
    /// @param broad_phase The broad phase algorithm to use.
    /// @param narrow_phase_ccd The narrow phase CCD algorithm to use.
    double compute_cfl_stepsize(
        const CollisionMesh& mesh,
        const std::vector<Pose>& poses_t0,
        const std::vector<Pose>& poses_t1,
        const double dhat,
        const double min_distance = 0.0,
        const std::shared_ptr<BroadPhase> broad_phase =
            make_default_broad_phase(),
        const NonlinearCCD& nonlinear_ccd = NonlinearCCD()) const;
};

} // namespace ipc::rigid
