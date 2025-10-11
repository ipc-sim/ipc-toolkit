#pragma once
#include <ipc/collision_mesh.hpp>
#include <ipc/smooth_contact/common.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {
class Primitive {
public:
    Primitive(const long id, const ParameterType& _param)
        : m_id(id)
        , param(_param)
    {
    }

    virtual ~Primitive() = default;

    bool operator==(const Primitive& other) const { return id() == other.id(); }

    long id() const { return m_id; }
    bool is_active() const { return m_is_active; }
    virtual int n_vertices() const = 0;
    virtual int n_dofs() const = 0;
    const std::vector<long>& vertex_ids() const { return m_vertex_ids; }

protected:
    std::vector<long> m_vertex_ids;
    long m_id;
    const ParameterType param;

    bool m_is_active = true;
};
} // namespace ipc
