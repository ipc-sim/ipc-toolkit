#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/smooth_contact/common.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class Primitive {
public:
    Primitive(const index_t id, const SmoothContactParameters& _params)
        : params(_params)
        , m_id(id)
    {
    }

    virtual ~Primitive() = default;

    bool operator==(const Primitive& other) const { return id() == other.id(); }

    index_t id() const { return m_id; }
    bool is_active() const { return m_is_active; }
    virtual int n_vertices() const = 0;
    virtual int n_dofs() const = 0;
    const std::vector<index_t>& vertex_ids() const { return m_vertex_ids; }

protected:
    const SmoothContactParameters params;
    std::vector<index_t> m_vertex_ids;
    index_t m_id;
    bool m_is_active = true;
};

} // namespace ipc
