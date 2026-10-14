#pragma once

#include "collisions/esp_collision.hpp"
#include "collisions/esp_collision_dict.hpp"

#include <ipc/collision_mesh.hpp>
#include <ipc/candidates/candidates.hpp>
#include <ipc/collisions/normal/edge_edge.hpp>

#include <memory>

namespace ipc {
class ESPCollisions {
public:
    /// @brief Create an empty collision set; build() fills it.
    ESPCollisions();
    virtual ~ESPCollisions();

    /// @brief Initialize the set of collisions used to compute the barrier potential.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param params ESP parameters.
    /// @param broad_phase Broad-phase method to use (default if nullptr).
    void build(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const ESPParameters& params,
        BroadPhase* broad_phase = nullptr);

    /// @brief Initialize the set of collisions used to compute the barrier potential.
    /// @param _candidates Distance candidates from which the collision set is
    ///        built. In 3D, build them with all_types = true.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param params ESP parameters.
    void build(
        const Candidates& _candidates,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const ESPParameters& params);

    // ------------------------------------------------------------------------

    /// @brief Get the number of collisions.
    size_t size() const;

    /// @brief Get if the collision set are empty.
    bool empty() const;

    /// @brief Clear the collision set.
    void clear();

    /// @brief Compute minimum distance between all contact candidates
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @return Squared minimum distance
    double compute_minimum_distance(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices) const;

    /// @brief Convert contact pairs to string. Vertex pairs also list their
    ///        squared distance, potential and gradient norm at vertices.
    std::string to_string(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const ESPParameters& params) const;

    /// @brief Number of contact candidates
    int n_candidates() const { return m_candidates.size(); }

public:
    /// @brief Collision candidates
    Candidates m_candidates;

    /// @brief Total number of collision pairs counted across all quadrature build functions
    size_t num_quadrature_collision_pairs = 0;

    /// @brief Collision sets of the quadrature points, keyed by primitive.
    /// @note Defined in the internal header esp_collision_maps.hpp, so that
    /// this header does not pull in Abseil, a private dependency.
    struct Maps;
    /// @brief Get the collision sets of the quadrature points.
    Maps& maps() { return *m_maps; }
    /// @brief Get the collision sets of the quadrature points.
    const Maps& maps() const { return *m_maps; }

private:
    std::unique_ptr<Maps> m_maps;
};
} // namespace ipc