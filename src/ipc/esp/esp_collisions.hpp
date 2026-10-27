#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/candidates/candidates.hpp>
#include <ipc/collisions/normal/edge_edge.hpp>
#include <ipc/esp/collisions/esp_collision.hpp>
#include <ipc/esp/collisions/esp_collision_dict.hpp>

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
    /// @throws std::invalid_argument if params.quad_order < 1 in 2D.
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
    /// @throws std::invalid_argument if params.quad_order < 1 in 2D.
    /// @throws std::invalid_argument In 3D, if the candidates were not built
    ///         with all_types = true (see Candidates::has_all_types).
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

    /// @brief Check that params use the quadrature this set was built with.
    ///
    /// The quadrature collision sets hold one entry per quadrature point of
    /// the rule used in build(): the Gauss-Lobatto edge rule of order
    /// params.quad_order in 2D, or params.quad_rule() in 3D. ESPPotential and
    /// the ESP friction build evaluate those entries with their own
    /// parameters, so a different rule would pair the wrong points with the
    /// wrong weights or read past the stored entries. A set that was never
    /// built (or was cleared) holds no entries and accepts any parameters.
    /// @param params Parameters about to be used with this set.
    /// @throws std::invalid_argument if params' quadrature differs from the
    ///         one used in build().
    void check_quadrature(const ESPParameters& params) const;

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

    /// @brief Mesh dimension of the last build(), or 0 if never built.
    int m_dim = 0;
    /// @brief Gauss-Lobatto edge quadrature order of the last 2D build().
    int m_quad_order = 0;
    /// @brief Face quadrature rule of the last 3D build().
    FaceQuadRule m_face_quad_rule;
    /// @brief Auto lattice spacing and cap of the last 3D build().
    double m_lattice_spacing = 0;
    int m_max_lattice_subdivisions = 0;
};
} // namespace ipc