#include "esp_collision_dict.hpp"

#include <array>

namespace ipc {
template <PointType pType, int DIM>
void ESPCollisionDict<pType, DIM>::initialize(
    const std::vector<index_t>& primitive_ids,
    const std::vector<index_t>& primary_vertex_ids,
    const unordered_map<std::array<index_t, 3>, std::shared_ptr<ESPCollision>>&
        map)
{
    assert(primary_vertex_ids.size() <= m_primary_vertex_ids.size());
    for (int i = 0; i < primary_vertex_ids.size(); i++) {
        m_primary_vertex_ids[i] = primary_vertex_ids[i];
    }

    assert(m_primitive_ids.size() <= m_primary_vertex_ids.size());
    for (int i = 0; i < primitive_ids.size(); i++) {
        m_primitive_ids[i] = primitive_ids[i];
    }

    std::set<index_t> vids;
    for (const auto& [key, val] : map) {
        for (const index_t vid : val->vertex_ids()) {
            vids.insert(vid);
        }
    }

    // Erase virtual vertex id, which is the largest in all ids
    if (pType != PointType::VERTEX && !map.empty()) {
        auto iter = std::prev(vids.end());
        auto ptr = map.begin().value();
        vids.erase(iter);
    }

    // Insert primary ids
    for (index_t vi : primary_vertex_ids) {
        vids.insert(vi);
    }

    m_vertex_ids.assign(vids.begin(), vids.end());
    assert(std::is_sorted(m_vertex_ids.begin(), m_vertex_ids.end()));

    for (int i = 0; i < m_vertex_ids.size(); i++) {
        m_vertex_ids_inverse[m_vertex_ids[i]] = i;
    }

    // Cache primary local ids
    for (int i = 0; i < m_primary_vertex_ids.size(); i++) {
        if (m_primary_vertex_ids[i] < 0) {
            break;
        }
        m_primary_local_ids[i] = vertex_ids_inverse(m_primary_vertex_ids[i]);
    }

    // Cache dofs
    m_dofs.resize(m_vertex_ids.size() * dim);
    for (int i = 0; i < m_vertex_ids.size(); i++) {
        for (int d = 0; d < dim; d++) {
            m_dofs[i * dim + d] = m_vertex_ids[i] * dim + d;
        }
    }

    // Cache primary dofs
    m_primary_dofs.clear();
    m_primary_dofs.reserve(m_primary_vertex_ids.size() * dim);
    for (index_t i : m_primary_vertex_ids) {
        if (i < 0) {
            break;
        }
        for (index_t d = 0; d < dim; d++) {
            m_primary_dofs.push_back(i * dim + d);
        }
    }

    // Convert unordered_map to typed vectors
    for (const auto& [key, val] : map) {
        switch (val->type()) {
        case ESPCollisionType::VERTEX_VERTEX:
            if constexpr (DIM == 2) {
                auto ptr = std::dynamic_pointer_cast<
                    ESPCollisionTemplate<Vertex2, Vertex2>>(val);
                assert(ptr);
                vv_collisions.push_back(*ptr);
            } else {
                auto ptr = std::dynamic_pointer_cast<
                    ESPCollisionTemplate<Vertex3, Vertex3>>(val);
                assert(ptr);
                vv_collisions.push_back(*ptr);
            }
            break;
        case ESPCollisionType::EDGE_VERTEX:
            if constexpr (DIM == 2) {
                auto ptr = std::dynamic_pointer_cast<
                    ESPCollisionTemplate<Vertex2, Edge2P1>>(val);
                assert(ptr);
                ev_collisions.push_back(*ptr);
            } else {
                auto ptr = std::dynamic_pointer_cast<
                    ESPCollisionTemplate<Edge3P1, Vertex3>>(val);
                assert(ptr);
                ev_collisions.push_back(*ptr);
            }
            break;
        case ESPCollisionType::FACE_VERTEX:
            if constexpr (DIM == 3) {
                auto ptr = std::dynamic_pointer_cast<
                    ESPCollisionTemplate<Face3P1, Vertex3>>(val);
                assert(ptr);
                fv_collisions.push_back(*ptr);
            } else {
                log_and_throw_error(
                    "FACE_VERTEX collision type not supported in 2D dict");
            }
            break;
        default:
            log_and_throw_error("Invalid collision type!");
        }
    }
}

template <PointType pType, int DIM>
ESPCollision& ESPCollisionDict<pType, DIM>::operator[](int i)
{
    return const_cast<ESPCollision&>(
        static_cast<const ESPCollisionDict&>(*this)[i]);
}

template <PointType pType, int DIM>
const ESPCollision& ESPCollisionDict<pType, DIM>::operator[](int i) const
{
    if (i < vv_collisions.size()) {
        return vv_collisions[i];
    } else {
        i -= vv_collisions.size();
        if (i < ev_collisions.size()) {
            return ev_collisions[i];
        } else {
            i -= ev_collisions.size();
            if (i < fv_collisions.size()) {
                return fv_collisions[i];
            } else {
                log_and_throw_error("Invalid index!");
            }
        }
    }
}

template <PointType pType, int DIM>
const std::vector<index_t>& ESPCollisionDict<pType, DIM>::vertex_ids() const
{
    return m_vertex_ids;
}

template <PointType pType, int DIM>
const std::vector<index_t>& ESPCollisionDict<pType, DIM>::primary_dofs() const
{
    return m_primary_dofs;
}

template <PointType pType, int DIM>
const std::vector<index_t>& ESPCollisionDict<pType, DIM>::dofs() const
{
    return m_dofs;
}

template <PointType pType, int DIM>
index_t ESPCollisionDict<pType, DIM>::vertex_ids_inverse(index_t id) const
{
    auto iter = m_vertex_ids_inverse.find(id);
    if (iter == m_vertex_ids_inverse.end()) {
        return -1;
    }
    return iter->second;
}

template class ESPCollisionDict<PointType::VERTEX>;
template class ESPCollisionDict<PointType::EDGE>;
template class ESPCollisionDict<PointType::FACE>;
template class ESPCollisionDict<PointType::EDGE, 2>;
template class ESPCollisionDict<PointType::VERTEX, 2>;
} // namespace ipc
