#pragma once

#include <ipc/collisions/collisions.hpp>
#include "vertex_vertex.hpp"
#include "edge_vertex.hpp"
#include "edge_edge.hpp"
#include "edge_edge_3d.hpp"
#include "face_vertex.hpp"
#include "face_face.hpp"

namespace ipc {

template <int dim>
class SmoothCollisions : public VirtualCollisions<(dim == 2) ? max_vert_2d : max_vert_3d> {
public:
    constexpr static int max_vert = (dim == 2) ? max_vert_2d : max_vert_3d;
    using Super = VirtualCollisions<max_vert>;
    /// @brief The type of the collisions.
    using value_type = SmoothCollision<max_vert>;

public:
    // SmoothCollisions() = default;
    SmoothCollisions(bool _use_high_order_quadrature = false)
    : use_high_order_quadrature(_use_high_order_quadrature)
    {
    }

    void compute_adaptive_dhat(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const ParameterType &param,
        const BroadPhaseMethod broad_phase_method = DEFAULT_BROAD_PHASE_METHOD);

    /// @brief Initialize the set of collisions used to compute the barrier potential.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param broad_phase_method Broad-phase method to use.
    void build(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const ParameterType &param,
        const bool use_adaptive_dhat,
        const BroadPhaseMethod broad_phase_method = DEFAULT_BROAD_PHASE_METHOD);

    /// @brief Initialize the set of collisions used to compute the barrier potential.
    /// @param candidates Distance candidates from which the collision set is built.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    void build(
        const Candidates& _candidates,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const ParameterType &param,
        const bool use_adaptive_dhat);

    // ------------------------------------------------------------------------

    /// @brief Get the number of collisions.
    size_t size() const override;

    /// @brief Get if the collision set are empty.
    bool empty() const override;

    /// @brief Clear the collision set.
    void clear() override;

    /// @brief Get a reference to collision at index i.
    /// @param i The index of the collision.
    /// @return A reference to the collision.
    value_type& operator[](size_t i) override;

    /// @brief Get a const reference to collision at index i.
    /// @param i The index of the collision.
    /// @return A const reference to the collision.
    const value_type& operator[](size_t i) const override;

    double compute_minimum_distance(
        const CollisionMesh& mesh, const Eigen::MatrixXd& vertices) const override;

    std::string
    to_string(const CollisionMesh& mesh, const Eigen::MatrixXd& vertices) const;

    void set_use_convergent_formulation(
        const bool use_convergent_formulation) override
    {
        logger().error("Smooth contact formulation doesn't have convergent version!");
    }

    void set_are_shape_derivatives_enabled(const bool are_shape_derivatives_enabled) override
    {
        logger().error("Smooth contact formulation doesn't have shape derivatives implemented!");
    }

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
        double out = std::max(vert_adaptive_dhat.maxCoeff(), edge_adaptive_dhat.maxCoeff());
        if constexpr (dim == 3)
            return std::max(out, face_adaptive_dhat.maxCoeff());
        return out;
    }

public:
    std::vector<std::shared_ptr<value_type>> collisions;

    const bool use_high_order_quadrature;

    Eigen::VectorXd vert_adaptive_dhat;
    Eigen::VectorXd edge_adaptive_dhat;
    Eigen::VectorXd face_adaptive_dhat;

    Candidates candidates;
};

} // namespace ipc
