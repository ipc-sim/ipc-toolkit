#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/candidates/candidates.hpp>
#include <ipc/collisions/normal/edge_edge.hpp>
#include <ipc/collisions/normal/edge_vertex.hpp>
#include <ipc/collisions/normal/face_vertex.hpp>
#include <ipc/collisions/normal/normal_collision.hpp>
#include <ipc/collisions/normal/plane_vertex.hpp>
#include <ipc/collisions/normal/vertex_vertex.hpp>
#include <ipc/smooth_contact/collisions/smooth_collision.hpp>

#include <Eigen/Core>

#include <vector>

namespace ipc {
class SmoothCollisions {
public:
    /// @brief The type of the collisions.
    using value_type = SmoothCollision;

public:
    SmoothCollisions() = default;
    virtual ~SmoothCollisions() = default;

    void compute_adaptive_dhat(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const SmoothContactParameters params,
        BroadPhase* broad_phase = nullptr);

    /// @brief Initialize the set of collisions used to compute the barrier potential.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param broad_phase_method Broad-phase method to use.
    void build(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const SmoothContactParameters params,
        const bool use_adaptive_dhat = false,
        BroadPhase* broad_phase = nullptr);

    /// @brief Initialize the set of collisions used to compute the barrier potential.
    /// @param candidates Distance candidates from which the collision set is built.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    void build(
        const Candidates& _candidates,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const SmoothContactParameters params,
        const bool use_adaptive_dhat = false);

    // ------------------------------------------------------------------------

    /// @brief Get the number of collisions.
    size_t size() const;

    /// @brief Get if the collision set are empty.
    bool empty() const;

    /// @brief Clear the collision set.
    void clear();

    /// @brief Get a reference to collision at index i.
    /// @param i The index of the collision.
    /// @return A reference to the collision.
    SmoothCollision& operator[](size_t i);

    /// @brief Get a const reference to collision at index i.
    /// @param i The index of the collision.
    /// @return A const reference to the collision.
    const SmoothCollision& operator[](size_t i) const;

    /// @brief Compute minimum distance between all contact candidates
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @return Squared minimum distance
    double compute_minimum_distance(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices) const;

    /// @brief Compute minimum distance between contact pairs with non-zero potential
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @return Squared minimum distance
    double compute_active_minimum_distance(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices) const;

    /// @brief Convert contact pairs to string
    std::string to_string(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const SmoothContactParameters& params) const;

    /// @brief Get per-vertex dhat value when dhat is adaptive
    double get_vert_dhat(int vert_id) const
    {
        if (vert_adaptive_dhat.size() > 1) {
            return vert_adaptive_dhat(vert_id);
        } else {
            return vert_adaptive_dhat(0);
        }
    }
    /// @brief Get per-edge dhat value when dhat is adaptive
    double get_edge_dhat(int edge_id) const
    {
        if (edge_adaptive_dhat.size() > 1) {
            return edge_adaptive_dhat(edge_id);
        } else {
            return edge_adaptive_dhat(0);
        }
    }
    /// @brief Get per-face dhat value when dhat is adaptive
    double get_face_dhat(int face_id) const
    {
        if (face_adaptive_dhat.size() > 1) {
            return face_adaptive_dhat(face_id);
        } else {
            return face_adaptive_dhat(0);
        }
    }
    /// @brief Get maximum dhat value when dhat is adaptive
    double get_max_dhat() const
    {
        double out = std::max(
            vert_adaptive_dhat.maxCoeff(), edge_adaptive_dhat.maxCoeff());
        if (face_adaptive_dhat.size() > 0) {
            return std::max(out, face_adaptive_dhat.maxCoeff());
        }
        return out;
    }

    /// @brief Number of contact candidates
    int n_candidates() const { return m_candidates.size(); }

public:
    /// @brief (active) collision pairs
    std::vector<std::shared_ptr<SmoothCollision>> collisions;

    /// @brief per-vertex adaptive dhat
    Eigen::VectorXd vert_adaptive_dhat;
    /// @brief per-edge adaptive dhat
    Eigen::VectorXd edge_adaptive_dhat;
    /// @brief per-face adaptive dhat
    Eigen::VectorXd face_adaptive_dhat;

    /// @brief Collision candidates
    Candidates m_candidates;
};

} // namespace ipc
