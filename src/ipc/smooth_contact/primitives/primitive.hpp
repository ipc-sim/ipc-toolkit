#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/smooth_contact/common.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/*
 * @brief Primitive class to be used in GCP contact pairs.
 *
 * Any derived class of Primitive requires not only the vertex positions on this
 * primitive, but also the closest point direction from this primitive to the
 * other. The constructor function will check if the primitive is active, and
 * the contact pair will be active only if both primitives are active.
 */
class Primitive {
public:
    Primitive(const index_t id, const SmoothContactParameters& _params)
        : params(_params)
        , m_id(id)
    {
    }

    virtual ~Primitive() = default;

    bool operator==(const Primitive& other) const { return id() == other.id(); }

    /// @brief Vertex/Edge/Face ID of this primitive
    index_t id() const { return m_id; }
    /// @brief Whether this primitive is active
    bool is_active() const { return m_is_active; }
    /// @brief Number of vertices on this primitive
    virtual int n_vertices() const = 0;
    /// @brief Number of vertices times the dimension
    virtual int n_dofs() const = 0;
    /// @brief Vertex IDs on this primitive
    const std::vector<index_t>& vertex_ids() const { return m_vertex_ids; }

protected:
    /// @brief GCP parameters
    const SmoothContactParameters params;
    /// @brief Vertex IDs on this primitive
    std::vector<index_t> m_vertex_ids;
    /// @brief Vertex/Edge/Face ID of this primitive
    index_t m_id;
    /// @brief Whether this primitive is active
    bool m_is_active = true;
};

} // namespace ipc
