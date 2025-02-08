#pragma once

#include <ipc/broad_phase/default_broad_phase.hpp>
#include <ipc/candidates/edge_edge.hpp>
#include <ipc/candidates/edge_vertex.hpp>
#include <ipc/candidates/face_vertex.hpp>
#include <ipc/candidates/vertex_vertex.hpp>

#include <Eigen/Core>

#include <vector>

namespace ipc {

class Candidates {
public:
    Candidates() = default;

    /// @brief Initialize the set of discrete collision detection candidates.
    /// @param mesh The surface of the collision mesh.
    /// @param vertices Surface vertex positions (rowwise).
    /// @param inflation_radius Amount to inflate the bounding boxes.
    /// @param broad_phase_method Broad phase method to use.
    void build(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const double inflation_radius = 0,
        const std::shared_ptr<BroadPhase> broad_phase =
            make_default_broad_phase());

    /// @brief Initialize the set of continuous collision detection candidates.
    /// @note Assumes the trajectory is linear.
    /// @param mesh The surface of the collision mesh.
    /// @param vertices_t0 Surface vertex starting positions (rowwise).
    /// @param vertices_t1 Surface vertex ending positions (rowwise).
    /// @param inflation_radius Amount to inflate the bounding boxes.
    /// @param broad_phase_method Broad phase method to use.
    void build(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const double inflation_radius = 0,
        const std::shared_ptr<BroadPhase> broad_phase =
            make_default_broad_phase());

    size_t size() const;

    bool empty() const;

    void clear();

    ContinuousCollisionCandidate& operator[](size_t i);
    const ContinuousCollisionCandidate& operator[](size_t i) const;

    /// @brief Determine if the step is collision free from the set of candidates.
    /// @note Assumes the trajectory is linear.
    /// @param mesh The collision mesh.
    /// @param vertices_t0 Surface vertex starting positions (rowwise).
    /// @param vertices_t1 Surface vertex ending positions (rowwise).
    /// @param min_distance The minimum distance allowable between any two elements.
    /// @param narrow_phase_ccd The narrow phase CCD algorithm to use.
    /// @returns True if <b>any</b> collisions occur.
    bool is_step_collision_free(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const double min_distance = 0.0,
        const NarrowPhaseCCD& narrow_phase_ccd =
            DEFAULT_NARROW_PHASE_CCD) const;

    /// @brief Computes a maximal step size that is collision free using the set of collision candidates.
    /// @note Assumes the trajectory is linear.
    /// @param mesh The collision mesh.
    /// @param vertices_t0 Surface vertex starting positions (rowwise). Assumed to be intersection free.
    /// @param vertices_t1 Surface vertex ending positions (rowwise).
    /// @param min_distance The minimum distance allowable between any two elements.
    /// @param narrow_phase_ccd The narrow phase CCD algorithm to use.
    /// @returns A step-size \f$\in [0, 1]\f$ that is collision free. A value of 1.0 if a full step and 0.0 is no step.
    double compute_collision_free_stepsize(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const double min_distance = 0.0,
        const NarrowPhaseCCD& narrow_phase_ccd =
            DEFAULT_NARROW_PHASE_CCD) const;

    /// @brief Computes a conservative bound on the largest-feasible step size for surface primitives not in collision.
    /// @param mesh The collision mesh.
    /// @param displacements Surface vertex displacements (rowwise).
    /// @param dhat Barrier activation distance.
    double compute_noncandidate_conservative_stepsize(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& displacements,
        const double dhat) const;

    /// @brief Computes a CFL-inspired CCD maximum step step size.
    /// @param mesh The collision mesh.
    /// @param vertices_t0 Surface vertex starting positions (rowwise).
    /// @param vertices_t1 Surface vertex ending positions (rowwise).
    /// @param dhat Barrier activation distance.
    /// @param min_distance The minimum distance allowable between any two elements.
    /// @param narrow_phase_ccd The narrow phase CCD algorithm to use.
    double compute_cfl_stepsize(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const double dhat,
        const double min_distance = 0.0,
        const std::shared_ptr<BroadPhase> broad_phase =
            make_default_broad_phase(),
        const NarrowPhaseCCD& narrow_phase_ccd =
            DEFAULT_NARROW_PHASE_CCD) const;

    bool save_obj(
        const std::string& filename,
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const;

public:
    std::vector<VertexVertexCandidate> vv_candidates;
    std::vector<EdgeVertexCandidate> ev_candidates;
    std::vector<EdgeEdgeCandidate> ee_candidates;
    std::vector<FaceVertexCandidate> fv_candidates;
};

} // namespace ipc
