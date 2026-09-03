#pragma once

#include <array>
#include <ipc/collision_mesh.hpp>
#include <ipc/esp/esp_parameters.hpp>
#include <ipc/math/span.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/**
 * @brief Base class for primitives used in ESP contact models.
 *
 * This class defines the common interface for geometric primitives (like
 * vertices and edges) involved in a ESP contact. Derived classes
 * are responsible for implementing the specific logic for their geometry type.
 */
class EspPrimitive {
public:
    constexpr static int MAX_NUM_VERTS = 3;
    EspPrimitive(const index_t id) : m_id(id) { }

    virtual ~EspPrimitive() = default;

    bool operator==(const EspPrimitive& other) const
    {
        return id() == other.id();
    }

    /// @brief Get the ID of this primitive (e.g., vertex ID, edge ID).
    index_t id() const { return m_id; }

    /// @brief Get the number of vertices in the primitive's stencil.
    virtual int n_vertices() const = 0;

    /// @brief Get the number of degrees of freedom for this primitive.
    virtual int n_dofs() const = 0;

    /// @brief Get the vertex IDs of the primitive's stencil.
    span<const index_t> vertex_ids() const
    {
        assert(MAX_NUM_VERTS >= n_vertices());
        return span<const index_t>(m_vertex_ids.data(), n_vertices());
    }

protected:
    /// @brief Vertex IDs of the stencil for this primitive.
    std::array<index_t, MAX_NUM_VERTS> m_vertex_ids { { -1, -1, -1 } };
    /// @brief The ID of this primitive.
    index_t m_id;
};

namespace {
    // Helper function to find the vertices adjacent to a given vertex in a 2D
    // mesh.
    std::vector<index_t>
    find_vertex_neighbors_2D(const CollisionMesh& mesh, const index_t v_id)
    {
        assert(mesh.dim() == 2);
        std::array<index_t, 2> neighbors;
        std::fill(neighbors.begin(), neighbors.end(), v_id);
        for (const auto& edge_id : mesh.vertex_edge_adjacencies()[v_id]) {
            const auto& edge = mesh.edges().row(edge_id);
            if (edge[0] == v_id) {
                neighbors[0] = edge[1];
            } else {
                neighbors[1] = edge[0];
            }
        }
        std::vector<index_t> neighbors_ordered;
        for (index_t n : neighbors) {
            if (n != v_id)
                neighbors_ordered.push_back(n);
        }
        return neighbors_ordered;
    }
} // namespace

/// @brief 2D vertex primitive with neighbor storage, for OGC.
class Vertex2ogc : public EspPrimitive {
public:
    static constexpr int N_CORE_POINTS = 1;
    static constexpr int N_POINTS = 1;
    static constexpr int DIM = 2;
    static constexpr int N_DOFS = N_POINTS * DIM;

    Vertex2ogc(
        const index_t id, const CollisionMesh& mesh, const Eigen::MatrixXd& V)
        : EspPrimitive(id)
    {
        n_verts = 0;
        m_vertex_ids[n_verts++] = id;
        std::vector<index_t> neighbors = find_vertex_neighbors_2D(mesh, id);
        for (const auto& neighbor_id : neighbors) {
            m_vertex_ids[n_verts++] = neighbor_id;
        }
    }

    int n_vertices() const override { return n_verts; }
    int n_dofs() const override { return n_vertices() * DIM; }

private:
    int n_verts;
};

class Edge2P1 : public EspPrimitive {
public:
    static constexpr int N_CORE_POINTS = 2;
    static constexpr int N_POINTS = 2;
    static constexpr int DIM = 2;
    static constexpr int N_DOFS = N_POINTS * DIM;

    Edge2P1(const index_t id, const CollisionMesh& mesh)
        : EspPrimitive(id)
    {
        m_vertex_ids[0] = mesh.edges()(id, 0);
        m_vertex_ids[1] = mesh.edges()(id, 1);
    }

    int n_vertices() const override { return 2; }
    int n_dofs() const override { return n_vertices() * DIM; }
};

/// @brief Simple 2D vertex primitive (single vertex, no neighbor storage).
class Vertex2 : public EspPrimitive {
public:
    static constexpr int N_CORE_POINTS = 1;
    static constexpr int N_POINTS = 1;
    static constexpr int DIM = 2;
    static constexpr int N_DOFS = N_POINTS * DIM;

    Vertex2(const index_t id, const CollisionMesh& /*mesh*/)
        : EspPrimitive(id)
    {
        m_vertex_ids[0] = id;
    }

    int n_vertices() const override { return 1; }
    int n_dofs() const override { return N_DOFS; }
};

class Vertex3 : public EspPrimitive {
public:
    static constexpr int N_CORE_POINTS = 1;
    static constexpr int N_POINTS = 1;
    static constexpr int DIM = 3;
    static constexpr int N_DOFS = N_POINTS * DIM;

    Vertex3(const index_t id, const CollisionMesh& mesh)
        : EspPrimitive(id)
    {
        m_vertex_ids[0] = id;
    }

    int n_vertices() const override { return 1; }
    int n_dofs() const override { return n_vertices() * DIM; }
};

class Edge3P1 : public EspPrimitive {
public:
    static constexpr int N_CORE_POINTS = 2;
    static constexpr int N_POINTS = 2;
    static constexpr int DIM = 3;
    static constexpr int N_DOFS = N_POINTS * DIM;

    Edge3P1(const index_t id, const CollisionMesh& mesh)
        : EspPrimitive(id)
    {
        m_vertex_ids[0] = mesh.edges()(id, 0);
        m_vertex_ids[1] = mesh.edges()(id, 1);
    }

    int n_vertices() const override { return 2; }
    int n_dofs() const override { return n_vertices() * DIM; }
};

class Face3P1 : public EspPrimitive {
public:
    static constexpr int N_CORE_POINTS = 3;
    static constexpr int N_POINTS = 3;
    static constexpr int DIM = 3;
    static constexpr int N_DOFS = N_POINTS * DIM;

    Face3P1(const index_t id, const CollisionMesh& mesh)
        : EspPrimitive(id)
    {
        m_vertex_ids[0] = mesh.faces()(id, 0);
        m_vertex_ids[1] = mesh.faces()(id, 1);
        m_vertex_ids[2] = mesh.faces()(id, 2);
    }

    int n_vertices() const override { return 3; }
    int n_dofs() const override { return n_vertices() * DIM; }
};

} // namespace ipc