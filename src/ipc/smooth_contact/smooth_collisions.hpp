#pragma once

#include <ipc/collisions/collisions.hpp>
#include "edge_vertex.hpp"
#include "edge_edge.hpp"
#include "face_vertex.hpp"
#include "face_face.hpp"

namespace ipc {

template <int dim, class TCollision>
class SmoothCollisions : public VirtualCollisions<2*dim> {
public:
    constexpr static int max_vert = 2 * dim;
    /// @brief The type of the collisions.
    using value_type = TCollision;

public:
    SmoothCollisions() = default;
    SmoothCollisions(bool _use_adaptive_eps)
    : use_adaptive_eps(_use_adaptive_eps)
    {
    }

    /// @brief Initialize the set of collisions used to compute the barrier potential.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param dhat The activation distance of the barrier.
    /// @param dmin Minimum distance.
    /// @param broad_phase_method Broad-phase method to use.
    void build(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const double dhat,
        const double dmin = 0,
        const BroadPhaseMethod broad_phase_method = DEFAULT_BROAD_PHASE_METHOD) override;

    /// @brief Initialize the set of collisions used to compute the barrier potential.
    /// @param candidates Distance candidates from which the collision set is built.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param dhat The activation distance of the barrier.
    /// @param  dmin  Minimum distance.
    void build(
        const Candidates& candidates,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const double dhat,
        const double dmin = 0) override;

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

public:
    std::vector<std::shared_ptr<value_type>> collisions;

    const bool use_adaptive_eps = false;
};

typedef SmoothCollisions<2, SmoothEdgeEdgeCollision<2>> SmoothCollisions2;
typedef SmoothCollisions<3, SmoothFaceFaceCollision> SmoothCollisions3;

} // namespace ipc
