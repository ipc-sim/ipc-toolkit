#pragma once

#include <ipc/smooth_contact/collisions/smooth_collision.hpp>

#include <ipc/collision_mesh.hpp>
#include <ipc/candidates/candidates.hpp>
#include <ipc/collisions/normal/edge_edge.hpp>
#include <ipc/collisions/normal/edge_vertex.hpp>
#include <ipc/collisions/normal/face_vertex.hpp>
#include <ipc/collisions/normal/normal_collision.hpp>
#include <ipc/collisions/normal/plane_vertex.hpp>
#include <ipc/collisions/normal/vertex_vertex.hpp>

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
        const ParameterType param,
        const std::shared_ptr<BroadPhase> broad_phase =
            make_default_broad_phase());

    /// @brief Initialize the set of collisions used to compute the barrier potential.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param broad_phase_method Broad-phase method to use.
    void build(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const ParameterType param,
        const bool use_adaptive_dhat = false,
        const std::shared_ptr<BroadPhase> broad_phase =
            make_default_broad_phase());

    /// @brief Initialize the set of collisions used to compute the barrier potential.
    /// @param candidates Distance candidates from which the collision set is built.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    void build(
        const Candidates& _candidates,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const ParameterType param,
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
    value_type& operator[](size_t i);

    /// @brief Get a const reference to collision at index i.
    /// @param i The index of the collision.
    /// @return A const reference to the collision.
    const value_type& operator[](size_t i) const;

    double compute_minimum_distance(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices) const;

    double compute_active_minimum_distance(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices) const;

    std::string to_string(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const ParameterType& params) const;

    double get_vert_dhat(int vert_id) const
    {
        if (vert_adaptive_dhat.size() > 1)
            return vert_adaptive_dhat(vert_id);
        else
            return vert_adaptive_dhat(0);
    }
    double get_edge_dhat(int edge_id) const
    {
        if (edge_adaptive_dhat.size() > 1)
            return edge_adaptive_dhat(edge_id);
        else
            return edge_adaptive_dhat(0);
    }
    double get_face_dhat(int face_id) const
    {
        if (face_adaptive_dhat.size() > 1)
            return face_adaptive_dhat(face_id);
        else
            return face_adaptive_dhat(0);
    }
    double get_max_dhat() const
    {
        double out = std::max(
            vert_adaptive_dhat.maxCoeff(), edge_adaptive_dhat.maxCoeff());
        if (face_adaptive_dhat.size() > 0)
            return std::max(out, face_adaptive_dhat.maxCoeff());
        return out;
    }

    inline int n_candidates() const { return candidates.size(); }

public:
    std::vector<std::shared_ptr<value_type>> collisions;

    Eigen::VectorXd vert_adaptive_dhat;
    Eigen::VectorXd edge_adaptive_dhat;
    Eigen::VectorXd face_adaptive_dhat;

    Candidates candidates;
};

} // namespace ipc
