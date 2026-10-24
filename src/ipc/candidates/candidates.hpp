#pragma once

#include <ipc/broad_phase/default_broad_phase.hpp>
#include <ipc/candidates/edge_edge.hpp>
#include <ipc/candidates/edge_vertex.hpp>
#include <ipc/candidates/face_vertex.hpp>
#include <ipc/candidates/plane_vertex.hpp>
#include <ipc/candidates/vertex_vertex.hpp>

#include <Eigen/Core>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

namespace ipc {

// Forward declarations
class AdjacencyList;

/// @brief A class for storing and managing collision candidates.
class Candidates {
public:
    Candidates() = default;

    /// @brief Initialize the set of discrete collision detection candidates.
    /// @param mesh The surface of the collision mesh.
    /// @param vertices Surface vertex positions (rowwise).
    /// @param inflation_radius Amount to inflate the bounding boxes.
    /// @param broad_phase Broad phase method to use.
    /// @param all_types Also detect the candidate types only ESP needs (see
    ///                  BroadPhase::detect_collision_candidates). Recorded in
    ///                  has_all_types.
    void build(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const double inflation_radius = 0,
        BroadPhase* broad_phase = nullptr,
        const bool all_types = false);

    /// @brief Initialize the set of continuous collision detection candidates.
    /// @note Assumes the trajectory is linear.
    /// @note These candidates never include the ESP-only types, so
    ///       has_all_types is false afterwards.
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

    /// @brief Get if the collision at i is a vertex-vertex collision.
    /// @param i The index of the collision.
    /// @return If the collision at i is a vertex-vertex collision.
    bool is_vertex_vertex(size_t i) const;

    /// @brief Get if the collision at i is an edge-vertex collision.
    /// @param i The index of the collision.
    /// @return If the collision at i is an edge-vertex collision.
    bool is_edge_vertex(size_t i) const;

    /// @brief Get if the collision at i is an edge-edge collision.
    /// @param i The index of the collision.
    /// @return If the collision at i is an edge-edge collision.
    bool is_edge_edge(size_t i) const;

    /// @brief Get if the collision at i is a face-vertex collision.
    /// @param i The index of the collision.
    /// @return If the collision at i is a face-vertex collision.
    bool is_face_vertex(size_t i) const;

    /// @brief Get if the collision at i is a plane-vertex collision.
    /// @param i The index of the collision.
    /// @return If the collision at i is a plane-vertex collision.
    bool is_plane_vertex(size_t i) const;

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

    /// @brief Build per-primitive adjacency lists from the candidates.
    /// @note Required before calling any of the *_set() queries below, which
    ///       return empty lists until then.
    /// @note The ve_set() and vf_set() results view storage inside this
    ///       object, so calling this again or clear() invalidates them.
    void convert_candidates_to_sets();

    /// @brief Vertices near vertex id (2D: also the endpoints of nearby edges).
    AdjacencyList vv_set(const CollisionMesh& mesh, index_t id) const;
    /// @brief Edges near vertex id.
    AdjacencyList ve_set(index_t id) const;
    /// @brief Faces near vertex id.
    AdjacencyList vf_set(index_t id) const;

    /// @brief Vertices near edge id, including its endpoints.
    AdjacencyList ev_set(const CollisionMesh& mesh, index_t id) const;
    /// @brief Edges near or adjacent to edge id, excluding id.
    AdjacencyList ee_set(const CollisionMesh& mesh, index_t id) const;
    /// @brief Faces near or adjacent to edge id, excluding faces containing id.
    AdjacencyList ef_set(const CollisionMesh& mesh, index_t id) const;

    /// @brief Vertices near face id, including its corners.
    AdjacencyList fv_set(const CollisionMesh& mesh, index_t id) const;
    /// @brief Edges near or adjacent to face id.
    AdjacencyList fe_set(const CollisionMesh& mesh, index_t id) const;
    /// @brief Faces near or adjacent to face id, excluding id.
    AdjacencyList ff_set(const CollisionMesh& mesh, index_t id) const;

public:
    std::vector<VertexVertexCandidate> vv_candidates;
    std::vector<EdgeVertexCandidate> ev_candidates;
    std::vector<EdgeEdgeCandidate> ee_candidates;
    std::vector<FaceVertexCandidate> fv_candidates;
    std::vector<PlaneVertexCandidate> pv_candidates;

    std::vector<EdgeFaceCandidate> ef_candidates;
    std::vector<FaceFaceCandidate> ff_candidates;

    /// @brief Whether these candidates include the types only ESP needs.
    ///
    /// In 3D, all_types adds whole-mesh vertex-vertex, edge-vertex, edge-face,
    /// and face-face candidates. ESPCollisions needs them and NormalCollisions
    /// would double count them, so both check this flag in 3D and throw on a
    /// mismatch. build() and clear() set it; if you fill the candidate lists
    /// by hand, set it to match what you filled.
    bool has_all_types = false;

private:
    static bool default_is_active(double candidate) { return true; }

    /// @brief Adjacency lists built by convert_candidates_to_sets().
    ///
    /// Held behind a shared pointer so this public header does not expose the
    /// storage layout, and so copying a Candidates shares the lists instead of
    /// copying their arrays.
    struct AdjacencySets;
    std::shared_ptr<AdjacencySets> m_sets;
};

/// @brief A sorted list of unique primitive ids, as returned by the adjacency
/// queries of Candidates.
///
/// The list either owns its ids or views ids stored inside the Candidates it
/// came from. We use this instead of std::set because every caller only
/// iterates, and a contiguous sorted range gives the same ascending order
/// without allocating a tree node per id. The tradeoff is that a view is
/// invalidated by clear() or convert_candidates_to_sets() on its Candidates.
class AdjacencyList {
public:
    AdjacencyList() = default;

    /// @brief Take ownership of ids that are already sorted and unique.
    explicit AdjacencyList(std::vector<index_t> sorted_unique_ids)
        : m_owned(std::move(sorted_unique_ids))
    {
    }

    /// @brief View ids owned elsewhere; they must stay alive and unchanged for
    /// as long as the list is used.
    static AdjacencyList view(const index_t* data, const size_t size)
    {
        AdjacencyList list;
        list.m_view = data;
        list.m_view_size = size;
        list.m_is_view = true;
        return list;
    }

    const index_t* begin() const { return m_is_view ? m_view : m_owned.data(); }
    const index_t* end() const { return begin() + size(); }
    size_t size() const { return m_is_view ? m_view_size : m_owned.size(); }
    bool empty() const { return size() == 0; }

    /// @brief Whether id is in the list (a binary search, since it is sorted).
    bool contains(const index_t id) const
    {
        return std::binary_search(begin(), end(), id);
    }

private:
    std::vector<index_t> m_owned;
    const index_t* m_view = nullptr;
    size_t m_view_size = 0;
    bool m_is_view = false;
};

} // namespace ipc
