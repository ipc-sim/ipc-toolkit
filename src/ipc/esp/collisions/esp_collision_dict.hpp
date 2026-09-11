#pragma once
#include "esp_collision_template.hpp"

#include <ipc/distance/distance_type.hpp>
#include <ipc/utils/unordered_map_and_set.hpp>

#include <array>
#include <cstdint>

namespace ipc {

enum class PointType : std::uint8_t { VERTEX, EDGE, FACE };

/// @brief A collection of collision pairs, they can be (Vert, Vert), (Vert, Edge), or (Vert, Face)
/// The first entry of the pairs is always "Vert", which could be actually:
///     1. A real vertex
///     2. A point on an edge, as the closest point between a pair of edges
///     3. A point at the face center / edge quadrature point
/// In 2 and 3, the "Vert" is a virtual vertex that does not exist in the
/// CollisionMesh, the ID of a virtual vertex is always #n_verts, i.e.
/// immediately after all real vertices.
/// @tparam DIM Spatial dimension (2 or 3). Default is 3.
template <PointType pType, int DIM = 3> class ESPCollisionDict {
public:
    static constexpr int DIMENSION = DIM;

    // Collision pair types depend on dimension.
    using VVType = std::conditional_t<
        DIM == 2,
        ESPCollisionTemplate<Vertex2, Vertex2>,
        ESPCollisionTemplate<Vertex3, Vertex3>>;
    using EVType = std::conditional_t<
        DIM == 2,
        ESPCollisionTemplate<Vertex2, Edge2P1>,
        ESPCollisionTemplate<Edge3P1, Vertex3>>;

    ESPCollisionDict() = default;
    ~ESPCollisionDict() = default;

    void initialize(
        const std::vector<index_t>& primitive_ids,
        const std::vector<index_t>& primary_vertex_ids,
        const unordered_map<
            std::array<index_t, 3>,
            std::shared_ptr<ESPCollision>>& map);

    const std::array<index_t, 4>& primary_vertex_ids() const
    {
        return m_primary_vertex_ids;
    }
    /// @brief Local indices of primary vertices (i.e. vertex_ids_inverse(primary_vertex_ids()[i]))
    const std::array<index_t, 4>& primary_local_ids() const
    {
        return m_primary_local_ids;
    }
    int size() const
    {
        return vv_collisions.size() + ev_collisions.size()
            + fv_collisions.size();
    }

    ESPCollision& operator[](int i);
    const ESPCollision& operator[](int i) const;

    template <
        PointType T = pType,
        typename = std::enable_if_t<T != PointType::EDGE>>
    int primitive_id() const
    {
        return m_primitive_ids[0];
    }

    template <
        PointType T = pType,
        typename = std::enable_if_t<T == PointType::EDGE>>
    std::array<index_t, 2> primitive_ids() const
    {
        return m_primitive_ids;
    }

    template <
        PointType T = pType,
        typename = std::enable_if_t<T == PointType::EDGE>>
    EdgeEdgeDistanceType ee_dtype() const
    {
        return m_ee_dtype;
    }

    template <
        PointType T = pType,
        typename = std::enable_if_t<T == PointType::EDGE>>
    void set_ee_dtype(EdgeEdgeDistanceType dtype)
    {
        m_ee_dtype = dtype;
    }

    /* These functions are only available after calling finish_insertion() */

    // Global indices of DoFs
    const std::vector<index_t>& dofs() const;
    const std::vector<index_t>& primary_dofs() const;
    // Global indices of vertices
    const std::vector<index_t>& vertex_ids() const;
    // Map from global vertex index to local vertex index
    index_t vertex_ids_inverse(index_t id) const;

private:
    std::vector<VVType> vv_collisions;
    std::vector<EVType> ev_collisions;
    std::vector<ESPCollisionTemplate<Face3P1, Vertex3>>
        fv_collisions; // unused in DIM=2

    std::array<index_t, 2> m_primitive_ids { { -1, -1 } };

    /// @brief Cached edge-edge distance type (only meaningful for PointType::EDGE)
    EdgeEdgeDistanceType m_ee_dtype = EdgeEdgeDistanceType::AUTO;

    /// @brief Primary vertices used to compute the virtual vertex
    /// When the quadrature point q is
    ///     - a vertex, this is that vertex id
    ///     - an edge point, this is the four vertices of the edge-edge pair
    ///     - a face point, this is the three vertices of the face
    /// When the size is smaller than 4, append -1 to entries not used.
    std::array<index_t, 4> m_primary_vertex_ids { { -1, -1, -1, -1 } };
    /// @brief Cached local indices: m_primary_local_ids[i] = vertex_ids_inverse(m_primary_vertex_ids[i])
    std::array<index_t, 4> m_primary_local_ids { { -1, -1, -1, -1 } };

    /// @brief Collection of all vertices in collision pairs, including the primary vertices, but not the virtual vertex
    std::vector<index_t> m_vertex_ids;
    /// @brief Inverse of m_vertex_ids
    std::map<index_t, index_t> m_vertex_ids_inverse;

    /// @brief Cached global DOF indices (computed once in initialize())
    std::vector<index_t> m_dofs;
    /// @brief Cached primary DOF indices (computed once in initialize())
    std::vector<index_t> m_primary_dofs;
};
} // namespace ipc