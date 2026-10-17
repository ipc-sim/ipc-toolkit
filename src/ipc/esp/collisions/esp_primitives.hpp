#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/esp/esp_parameters.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <array>

namespace ipc {

/**
 * @brief Base class for primitives used in ESP contact models.
 *
 * This class defines the common interface for geometric primitives (like
 * vertices and edges) involved in a ESP contact. Derived classes
 * are responsible for implementing the specific logic for their geometry type.
 */
class ESPPrimitive {
public:
    constexpr static int MAX_NUM_VERTS = 3;
    ESPPrimitive(const index_t id) : m_id(id) { }

    virtual ~ESPPrimitive() = default;

    bool operator==(const ESPPrimitive& other) const
    {
        return id() == other.id();
    }

    /// @brief Get the ID of this primitive (e.g., vertex ID, edge ID).
    index_t id() const { return m_id; }

    /// @brief Get the number of vertices in the primitive's stencil.
    virtual int n_vertices() const = 0;

    /// @brief Get the number of degrees of freedom for this primitive.
    virtual int n_dofs() const = 0;

    /// @brief Get the i-th vertex ID of the primitive's stencil.
    index_t vertex_id(const index_t i) const
    {
        assert(i >= 0 && i < n_vertices());
        return m_vertex_ids[i];
    }

protected:
    /// @brief Vertex IDs of the stencil for this primitive.
    std::array<index_t, MAX_NUM_VERTS> m_vertex_ids { { -1, -1, -1 } };
    /// @brief The ID of this primitive.
    index_t m_id;
};

class Edge2P1 : public ESPPrimitive {
public:
    static constexpr int N_CORE_POINTS = 2;
    static constexpr int N_POINTS = 2;
    static constexpr int DIM = 2;
    static constexpr int N_DOFS = N_POINTS * DIM;

    Edge2P1(const index_t id, const CollisionMesh& mesh) : ESPPrimitive(id)
    {
        m_vertex_ids[0] = mesh.edges()(id, 0);
        m_vertex_ids[1] = mesh.edges()(id, 1);
    }

    int n_vertices() const override { return 2; }
    int n_dofs() const override { return n_vertices() * DIM; }
};

/// @brief Simple 2D vertex primitive (single vertex, no neighbor storage).
class Vertex2 : public ESPPrimitive {
public:
    static constexpr int N_CORE_POINTS = 1;
    static constexpr int N_POINTS = 1;
    static constexpr int DIM = 2;
    static constexpr int N_DOFS = N_POINTS * DIM;

    Vertex2(const index_t id, const CollisionMesh& /*mesh*/) : ESPPrimitive(id)
    {
        m_vertex_ids[0] = id;
    }

    int n_vertices() const override { return 1; }
    int n_dofs() const override { return N_DOFS; }
};

class Vertex3 : public ESPPrimitive {
public:
    static constexpr int N_CORE_POINTS = 1;
    static constexpr int N_POINTS = 1;
    static constexpr int DIM = 3;
    static constexpr int N_DOFS = N_POINTS * DIM;

    Vertex3(const index_t id, const CollisionMesh& mesh) : ESPPrimitive(id)
    {
        m_vertex_ids[0] = id;
    }

    int n_vertices() const override { return 1; }
    int n_dofs() const override { return n_vertices() * DIM; }
};

class Edge3P1 : public ESPPrimitive {
public:
    static constexpr int N_CORE_POINTS = 2;
    static constexpr int N_POINTS = 2;
    static constexpr int DIM = 3;
    static constexpr int N_DOFS = N_POINTS * DIM;

    Edge3P1(const index_t id, const CollisionMesh& mesh) : ESPPrimitive(id)
    {
        m_vertex_ids[0] = mesh.edges()(id, 0);
        m_vertex_ids[1] = mesh.edges()(id, 1);
    }

    int n_vertices() const override { return 2; }
    int n_dofs() const override { return n_vertices() * DIM; }
};

class Face3P1 : public ESPPrimitive {
public:
    static constexpr int N_CORE_POINTS = 3;
    static constexpr int N_POINTS = 3;
    static constexpr int DIM = 3;
    static constexpr int N_DOFS = N_POINTS * DIM;

    Face3P1(const index_t id, const CollisionMesh& mesh) : ESPPrimitive(id)
    {
        m_vertex_ids[0] = mesh.faces()(id, 0);
        m_vertex_ids[1] = mesh.faces()(id, 1);
        m_vertex_ids[2] = mesh.faces()(id, 2);
    }

    int n_vertices() const override { return 3; }
    int n_dofs() const override { return n_vertices() * DIM; }
};

} // namespace ipc