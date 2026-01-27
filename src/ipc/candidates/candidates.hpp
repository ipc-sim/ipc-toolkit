#pragma once

#include <ipc/broad_phase/default_broad_phase.hpp>
#include <ipc/candidates/edge_edge.hpp>
#include <ipc/candidates/edge_vertex.hpp>
#include <ipc/candidates/face_vertex.hpp>
#include <ipc/candidates/vertex_vertex.hpp>

#include <Eigen/Core>

#include <vector>

namespace ipc {

/// @brief A class for storing and managing collision candidates.
class Candidates {
public:
    Candidates() = default;

    /// @brief Initialize the set of discrete collision detection candidates.
    /// @param mesh The surface of the collision mesh.
    /// @param vertices Surface vertex positions (rowwise).
    /// @param inflation_radius Amount to inflate the bounding boxes.
    /// @param broad_phase Broad phase method to use.
    void build(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const double inflation_radius = 0,
        BroadPhase* broad_phase = nullptr);

    /// @brief Initialize the set of continuous collision detection candidates.
    /// @note Assumes the trajectory is linear.
    /// @param mesh The surface of the collision mesh.
    /// @param vertices_t0 Surface vertex starting positions (rowwise).
    /// @param vertices_t1 Surface vertex ending positions (rowwise).
    /// @param inflation_radius Amount to inflate the bounding boxes.
    /// @param broad_phase Broad phase method to use.
    void build(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
        Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
        const double inflation_radius = 0,
        BroadPhase* broad_phase = nullptr);

    /// @brief Get the number of collision candidates.
    /// @return The number of collision candidates.
    size_t size() const;

    /// @brief Check if there are no collision candidates.
    /// @return True if there are no collision candidates, false otherwise.
    bool empty() const;

    /// @brief Clear all collision candidates.
    void clear();

    /// @brief Get a collision stencil by index.
    /// @param i The index of the collision stencil.
    /// @return A reference to the collision stencil.
    CollisionStencil& operator[](size_t i);

    /// @brief Get a collision stencil by index.
    /// @param i The index of the collision stencil.
    /// @return A const reference to the collision stencil.
    const CollisionStencil& operator[](size_t i) const;

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
        Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
        Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
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
        Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
        Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
        const double min_distance = 0.0,
        const NarrowPhaseCCD& narrow_phase_ccd =
            DEFAULT_NARROW_PHASE_CCD) const;

    /// @brief Computes a conservative bound on the largest-feasible step size for surface primitives not in collision.
    /// @param mesh The collision mesh.
    /// @param displacements Surface vertex displacements (rowwise).
    /// @param dhat Barrier activation distance.
    /// @return A step-size \f$\in [0, 1]\f$ that is collision free for non-candidate elements.
    double compute_noncandidate_conservative_stepsize(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> displacements,
        const double dhat) const;

    /// @brief Computes a CFL-inspired CCD maximum step step size.
    /// @param mesh The collision mesh.
    /// @param vertices_t0 Surface vertex starting positions (rowwise).
    /// @param vertices_t1 Surface vertex ending positions (rowwise).
    /// @param dhat Barrier activation distance.
    /// @param min_distance The minimum distance allowable between any two elements.
    /// @param broad_phase The broad phase algorithm to use.
    /// @param narrow_phase_ccd The narrow phase CCD algorithm to use.
    /// @returns A step-size \f$\in [0, 1]\f$ that is collision free.
    double compute_cfl_stepsize(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
        Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
        const double dhat,
        const double min_distance = 0.0,
        BroadPhase* broad_phase = nullptr,
        const NarrowPhaseCCD& narrow_phase_ccd =
            DEFAULT_NARROW_PHASE_CCD) const;

    /// @brief Compute the maximum distance every vertex can move (independently) without colliding with any other element.
    /// @note Cap the value at the inflation radius used to build the candidates.
    /// @param mesh The collision mesh.
    /// @param vertices Collision mesh vertex positions (rowwise).
    /// @param inflation_radius The inflation radius used to build the candidates.
    /// @param min_distance The minimum distance allowable between any two elements.
    /// @return A vector of minimum distances, one for each vertex.
    Eigen::VectorXd compute_per_vertex_safe_distances(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const double inflation_radius,
        const double min_distance = 0.0) const;

    // clang-format off
    // /// @brief Compute the maximum distance every vertex can move (independently) without colliding with any other element.
    // /// @note Cap the value at one.
    // /// @param mesh The collision mesh.
    // /// @param vertices_t0 Surface vertex starting positions (rowwise).
    // /// @param vertices_t1 Surface vertex ending positions (rowwise).
    // /// @param min_distance The minimum distance allowable between any two elements.
    // /// @return A vector of values in [0, 1], one for each vertex.
    // Eigen::VectorXd compute_per_vertex_collision_free_stepsize(
    //     const CollisionMesh& mesh,
    //     Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
    //     Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
    //     const double min_distance = 0.0,
    //     const NarrowPhaseCCD& narrow_phase_ccd =
    //         DEFAULT_NARROW_PHASE_CCD) const;
    // clang-format on

    // == Convert to subelement candidates ====================================

    /// @brief Convert edge-vertex candidates to vertex-vertex candidates.
    /// @param mesh The collision mesh.
    /// @param vertices Collision mesh vertex positions (rowwise).
    /// @param is_active (Optional) Function to determine if a candidate is active.
    /// @return Vertex-vertex candidates derived from edge-vertex candidates.
    std::vector<VertexVertexCandidate> edge_vertex_to_vertex_vertex(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const std::function<bool(double)>& is_active = default_is_active) const;

    /// @brief Convert face-vertex candidates to vertex-vertex candidates.
    /// @param mesh The collision mesh.
    /// @param vertices Collision mesh vertex positions (rowwise).
    /// @param is_active (Optional) Function to determine if a candidate is active.
    /// @return Vertex-vertex candidates derived from face-vertex candidates.
    std::vector<VertexVertexCandidate> face_vertex_to_vertex_vertex(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const std::function<bool(double)>& is_active = default_is_active) const;

    /// @brief Convert face-vertex candidates to edge-vertex candidates.
    /// @param mesh The collision mesh.
    /// @param vertices Collision mesh vertex positions (rowwise).
    /// @param is_active (Optional) Function to determine if a candidate is active.
    /// @return Edge-vertex candidates derived from face-vertex candidates.
    std::vector<EdgeVertexCandidate> face_vertex_to_edge_vertex(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const std::function<bool(double)>& is_active = default_is_active) const;

    /// @brief Convert edge-edge candidates to edge-vertex candidates.
    /// @param mesh The collision mesh.
    /// @param vertices Collision mesh vertex positions (rowwise).
    /// @param is_active (Optional) Function to determine if a candidate is active.
    /// @return Edge-vertex candidates derived from edge-edge candidates.
    std::vector<EdgeVertexCandidate> edge_edge_to_edge_vertex(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const std::function<bool(double)>& is_active = default_is_active) const;

    // == Save candidates to file =============================================

    /// @brief Write the collision candidates to an OBJ file.
    /// @param filename The name of the file to write the candidates to.
    /// @param vertices Collision mesh vertex positions (rowwise).
    /// @param edges Collision mesh edge indices (rowwise).
    /// @param faces Collision mesh face indices (rowwise).
    /// @return True if the file was written successfully.
    bool write_obj(
        const std::string& filename,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const;

public:
    std::vector<VertexVertexCandidate> vv_candidates;
    std::vector<EdgeVertexCandidate> ev_candidates;
    std::vector<EdgeEdgeCandidate> ee_candidates;
    std::vector<FaceVertexCandidate> fv_candidates;

private:
    static bool default_is_active(double candidate) { return true; }
};

} // namespace ipc
