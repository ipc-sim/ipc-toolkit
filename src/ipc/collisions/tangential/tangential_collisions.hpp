#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/collisions/tangential/edge_edge.hpp>
#include <ipc/collisions/tangential/edge_vertex.hpp>
#include <ipc/collisions/tangential/face_vertex.hpp>
#include <ipc/collisions/tangential/vertex_vertex.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <map>
#include <functional>

namespace ipc {

class TangentialCollisions {
public:
    /// @brief The type of the collisions.
    using value_type = TangentialCollision;

    /// @brief Blend type enum for selecting the method of blending friction coefficients.
    enum class BlendType {
        MAX,
        HARMONIC_MEAN
    };

    /// @brief The friction coefficients for a pair of materials.
    struct MaterialPairFriction {
        double s_mu;
        double k_mu;
    };

public:
    TangentialCollisions() = default;

    void build(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const NormalCollisions& collisions,
        const NormalPotential& normal_potential,
        double normal_stiffness,
        double mu)
    {
        this->build(
            mesh, vertices, collisions, normal_potential, normal_stiffness,
            Eigen::VectorXd::Constant(vertices.rows(), mu));
    }

    void build(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const NormalCollisions& collisions,
        const NormalPotential& normal_potential,
        const double normal_stiffness,
        Eigen::ConstRef<Eigen::VectorXd> mus,
        const std::function<double(double, double)>& blend_mu = blend_mu_wrapper);

    void build(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const NormalCollisions& collisions,
        const NormalPotential& normal_potential,
        double barrier_stiffness,
        double mu,
        double s_mu,
        double k_mu);

    // ------------------------------------------------------------------------

    /// @brief Get the number of friction collisions.
    size_t size() const;

    /// @brief Get if the friction collisions are empty.
    bool empty() const;

    /// @brief Clear the friction collisions.
    void clear();

    /// @brief Get a reference to collision at index i.
    TangentialCollision& operator[](const size_t i);

    /// @brief Get a const reference to collision at index i.
    const TangentialCollision& operator[](const size_t i) const;

    static double default_blend_mu(
        double mu1, double mu2, BlendType blend_type = BlendType::MAX)
    {
        if (blend_type == BlendType::MAX) {
            return std::max(mu1, mu2);
        } else {
            // HARMONIC_MEAN
            if (mu1 <= 0 || mu2 <= 0) {
                return 0;
            }
            return 2 * (mu1 * mu2) / (mu1 + mu2);
        }
    }

    static double blend_mu_wrapper(double mu1, double mu2) {
        return default_blend_mu(mu1, mu2);
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

    /// @brief The blend type to use when combining friction coefficients.
    BlendType blend_type = BlendType::MAX;
};

} // namespace ipc
