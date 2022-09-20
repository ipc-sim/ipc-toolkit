#pragma once

#include <ipc/collisions/plane_vertex.hpp>

#include <Eigen/Core>

#include <vector>

namespace ipc {

/// @brief Construct a set of point-plane distance constraints used to compute
/// the barrier potential.
///
/// @note The given pv_constraints will be cleared.
///
/// @param[in]  V  Vertex positions as rows of a matrix.
/// @param[in]  plane_origins  Plane origins as rows of a matrix.
/// @param[in]  plane_normals  Plane normals as rows of a matrix.
/// @param[in]  dhat  The activation distance of the barrier.
/// @param[out] pv_constraints  The constructed set of constraints.
/// @param[in]  dmin  Minimum distance.
/// @param[in] can_collide
///     A function that takes a vertex ID (row numbers in V) and a plane ID (row
///     number in plane_origins) then returns true if the vertex can collide
///     with the plane. By default all points can collide with all planes.
void construct_point_plane_constraint_set(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXd& plane_origins,
    const Eigen::MatrixXd& plane_normals,
    const double dhat,
    std::vector<PlaneVertexConstraint>& pv_constraints,
    const double dmin = 0,
    const std::function<bool(size_t, size_t)>& can_collide =
        [](size_t, size_t) { return true; });

///////////////////////////////////////////////////////////////////////////////
// Collision detection

/// @brief Determine if the step is collision free.
///
/// @note Assumes the trajectory is linear.
///
/// @param[in] V0 Vertex positions at start as rows of a matrix.
/// @param[in] V1 Vertex positions at end as rows of a matrix.
/// @param[in] plane_origins  Plane origins as rows of a matrix.
/// @param[in] plane_normals  Plane normals as rows of a matrix.
/// @param[in] can_collide
///     A function that takes a vertex ID (row numbers in V) and a plane ID (row
///     number in plane_origins) then returns true if the vertex can collide
///     with the plane. By default all points can collide with all planes.
/// @returns True if <b>any</b> collisions occur.
bool is_step_point_plane_collision_free(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXd& plane_origins,
    const Eigen::MatrixXd& plane_normals,
    const std::function<bool(size_t, size_t)>& can_collide =
        [](size_t, size_t) { return true; });

/// @brief Computes a maximal step size that is collision free.
///
/// @note Assumes V0 is intersection free.
/// @note Assumes the trajectory is linear.
/// @note A value of 1.0 if a full step and 0.0 is no step.
///
/// @param[in] V0 Vertex positions at start as rows of a matrix.
/// @param[in] V1 Vertex positions at end as rows of a matrix.
/// @param[in] plane_origins  Plane origins as rows of a matrix.
/// @param[in] plane_normals  Plane normals as rows of a matrix.
/// @param[in] can_collide
///     A function that takes a vertex ID (row numbers in V) and a plane ID (row
///     number in plane_origins) then returns true if the vertex can collide
///     with the plane. By default all points can collide with all planes.
/// @returns A step-size \f$\in [0, 1]\f$ that is collision free.
double compute_point_plane_collision_free_stepsize(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXd& plane_origins,
    const Eigen::MatrixXd& plane_normals,
    const std::function<bool(size_t, size_t)>& can_collide =
        [](size_t, size_t) { return true; });

} // namespace ipc
