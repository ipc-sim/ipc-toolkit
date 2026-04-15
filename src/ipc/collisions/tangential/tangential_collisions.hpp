#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/collisions/tangential/edge_edge.hpp>
#include <ipc/collisions/tangential/edge_vertex.hpp>
#include <ipc/collisions/tangential/face_vertex.hpp>
#include <ipc/collisions/tangential/plane_vertex.hpp>
#include <ipc/collisions/tangential/tangential_collision.hpp>
#include <ipc/collisions/tangential/vertex_vertex.hpp>
#include <ipc/smooth_contact/smooth_collisions.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ipc {

class TangentialCollisions {
public:
    /// @brief The type of the collisions.
    using value_type = TangentialCollision;

public:
    TangentialCollisions() = default;

    /// @brief Build the tangential collisions.
    /// @param mesh The collision mesh.
    /// @param vertices The vertices of the mesh.
    /// @param collisions The set of normal collisions.
    /// @param normal_potential The normal potential.
    /// @param mu The coefficient of friction.
    void build(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const NormalCollisions& collisions,
        const NormalPotential& normal_potential,
        double mu)
    {
        this->build(mesh, vertices, collisions, normal_potential, mu, mu);
    }

    /// @brief Build the tangential collisions.
    /// @param mesh The collision mesh.
    /// @param vertices The vertices of the mesh.
    /// @param collisions The set of normal collisions.
    /// @param normal_potential The normal potential.
    /// @param mu_s The static friction coefficient.
    /// @param mu_k The kinetic friction coefficient.
    void build(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const NormalCollisions& collisions,
        const NormalPotential& normal_potential,
        double mu_s,
        double mu_k)
    {
        this->build(
            mesh, vertices, collisions, normal_potential,
            Eigen::VectorXd::Constant(vertices.rows(), mu_s),
            Eigen::VectorXd::Constant(vertices.rows(), mu_k));
    }

    /// @brief Build the tangential collisions.
    /// @param mesh The collision mesh.
    /// @param vertices The vertices of the mesh.
    /// @param collisions The set of normal collisions.
    /// @param normal_potential The normal potential.
    /// @param mu_k The kinetic friction coefficient per vertex.
    /// @param mu_s The static friction coefficient per vertex.
    /// @param blend_mu Function to blend vertex-based coefficients of friction.
    void build(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const NormalCollisions& collisions,
        const NormalPotential& normal_potential,
        Eigen::ConstRef<Eigen::VectorXd> mu_s,
        Eigen::ConstRef<Eigen::VectorXd> mu_k,
        const std::function<double(double, double)>& blend_mu =
            default_blend_mu);

    /// @brief Build the tangential collisions for smooth contact.
    /// @param mesh The collision mesh.
    /// @param vertices The vertices of the mesh.
    /// @param collisions The set of normal collisions.
    /// @param params Parameters of Geometric Contact Potential
    /// @param normal_stiffness Stiffness of the normal potential.
    /// @param mu_s The static friction coefficient per vertex.
    /// @param mu_k The kinetic friction coefficient per vertex.
    /// @param blend_mu Function to blend vertex-based coefficients of friction. Defaults to average.
    void build(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const SmoothCollisions& collisions,
        const SmoothContactParameters& params,
        const double normal_stiffness,
        Eigen::ConstRef<Eigen::VectorXd> mu_s,
        Eigen::ConstRef<Eigen::VectorXd> mu_k,
        const std::function<double(double, double)>& blend_mu =
            default_blend_mu);

    /// @brief Set lagged effective μ to scalar mu_s/mu_k on every collision (after build).
    void reset_lagged_anisotropic_friction_coefficients();

    /// @brief Refresh lagged matchstick effective μ from lagged geometry and current slip.
    /// @param mesh Collision mesh (edges/faces for DOF gathering).
    /// @param rest_positions Rest configuration (rows = vertices).
    /// @param lagged_displacements Displacements at lagged state (same shape as rest).
    /// @param velocities World velocities at stencil vertices (same shape).
    /// @note Call after build and after any lagged-geometry or velocity update
    ///       used for slip (typically each Newton iteration). Friction
    ///       force/gradient/Hessian/Jacobians then use these cached scalars and
    ///       do not differentiate μ_eff with respect to slip in that
    ///       evaluation. Directional lagged μ is currently activated when
    ///       mu_s_aniso is nonzero in a 2D tangent space (3D simulation).
    void update_lagged_anisotropic_friction_coefficients(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> rest_positions,
        Eigen::ConstRef<Eigen::MatrixXd> lagged_displacements,
        Eigen::ConstRef<Eigen::MatrixXd> velocities);

    // ------------------------------------------------------------------------

    /// @brief Get the number of friction collisions.
    size_t size() const;

    /// @brief Get if the friction collisions are empty.
    bool empty() const;

    /// @brief Clear the friction collisions.
    void clear();

    /// @brief Get a reference to collision at index i.
    /// @param i The index of the collision.
    /// @return A reference to the collision.
    TangentialCollision& operator[](const size_t i);

    /// @brief Get a const reference to collision at index i.
    /// @param i The index of the collision.
    /// @return A const reference to the collision.
    const TangentialCollision& operator[](const size_t i) const;

    static double default_blend_mu(double mu0, double mu1)
    {
        // return mu0 * mu1;
        // return std::min(mu0, mu1);
        // return std::max(mu0, mu1);
        return (mu0 + mu1) / 2;
    }

public:
    /// @brief Vertex-vertex tangential collisions.
    std::vector<VertexVertexTangentialCollision> vv_collisions;
    /// @brief Edge-vertex tangential collisions.
    std::vector<EdgeVertexTangentialCollision> ev_collisions;
    /// @brief Edge-edge tangential collisions.
    std::vector<EdgeEdgeTangentialCollision> ee_collisions;
    /// @brief Face-vertex tangential collisions.
    std::vector<FaceVertexTangentialCollision> fv_collisions;
    /// @brief Plane-vertex tangential collisions.
    std::vector<PlaneVertexTangentialCollision> pv_collisions;
};

} // namespace ipc
