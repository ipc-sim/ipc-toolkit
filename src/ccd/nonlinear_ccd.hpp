#pragma once

#include <ipc/ccd/ccd.hpp>

#include <filib/interval.hpp>

#include <functional>

namespace ipc {

/// Find time-of-impact between two rigid bodies
// bool compute_piecewise_linear_edge_vertex_time_of_impact(
//     const RigidBody& bodyA,
//     const PoseD& poseA_t0, // Pose of bodyA at t=0
//     const PoseD& poseA_t1, // Pose of bodyA at t=1
//     size_t vertex_id,      // In bodyA
//     const RigidBody& bodyB,
//     const PoseD& poseB_t0, // Pose of bodyB at t=0
//     const PoseD& poseB_t1, // Pose of bodyB at t=1
//     size_t edge_id,        // In bodyB
//     double& toi,
//     double earliest_toi = 1, // Only search for collision in [0,
//     earliest_toi], double minimum_separation_distance = 0, double
//     toi_tolerance = Constants::RIGID_CCD_TOI_TOL);

bool edge_edge_nonlinear_ccd(
    const std::function<interval(const interval&)>& ea0, // ea0([t0, t1])
    const std::function<interval(const interval&)>& ea1, // ea1([t0, t1])
    const std::function<interval(const interval&)>& eb0, // eb0([t0, t1])
    const std::function<interval(const interval&)>& eb1, // eb1([t0, t1])
    double& toi,
    const double tmax = 1.0,
    const double min_distance = 0,
    const double tolerance = DEFAULT_CCD_TOLERANCE,
    const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
    const double conservative_rescaling = DEFAULT_CCD_CONSERVATIVE_RESCALING);

/// Find time-of-impact between two rigid bodies
// bool compute_piecewise_linear_face_vertex_time_of_impact(
//     const RigidBody& bodyA,
//     const PoseD& poseA_t0, // Pose of bodyA at t=0
//     const PoseD& poseA_t1, // Pose of bodyA at t=1
//     size_t vertex_id,      // In bodyA
//     const RigidBody& bodyB,
//     const PoseD& poseB_t0, // Pose of bodyB at t=0
//     const PoseD& poseB_t1, // Pose of bodyB at t=1
//     size_t face_id,        // In bodyB
//     double& toi,
//     double earliest_toi = 1, // Only search for collision in [0,
//     earliest_toi], double minimum_separation_distance = 0, double
//     toi_tolerance = Constants::RIGID_CCD_TOI_TOL);

} // namespace ipc
