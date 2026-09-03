#include "adaptive_support.hpp"

#include "collisions/esp_quadrature.hpp"
#include "collisions/vertex_matrix_view.hpp"
#include "esp_collisions.hpp"

#include <ipc/candidates/candidates.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/gcp/distance/edge_edge.hpp>
#include <ipc/utils/logger.hpp>

namespace ipc {

AdaptiveSupport::AdaptiveSupport(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> rest_positions,
    const EspParameters& params)
    : m_mesh(&mesh)
{
    const int nv = mesh.num_vertices();
    m_values.setConstant(nv, params.dhat);

    EspCollisions collisions;
    collisions.build(mesh, rest_positions, params);

    if (collisions.empty())
        return;

    // Returns the mesh vertex IDs in a collision pair that belong to the
    // PRIMITIVE (i.e., not the source quadrature point). Source vertices are
    // those listed in the dict's primary_vertex_ids. Virtual vertices
    // (id >= nv) are skipped.
    auto get_primitive_vids =
        [&](const EspCollision& cc) -> std::vector<index_t> {
        std::vector<index_t> pvids;
        for (int i = 0; i < cc.num_vertices(); i++) {
            const index_t vid = cc.vertex_id(i);
            if (vid < static_cast<index_t>(nv)) {
                pvids.push_back(vid);
            }
        }
        return pvids;
    };

    if (mesh.dim() == 3) {
        struct ActivePair {
            const EspCollision* cc;
            bool needs_extended;
            Eigen::RowVector3d qp_pos;
            std::vector<index_t> primitive_vids;
        };

        std::vector<ActivePair> active_pairs;
        const auto& face_quad_rule = params.get_quad_rule();

        for (index_t f = 0; f < static_cast<index_t>(mesh.num_faces()); f++) {
            // Face interior quadrature points
            if (!face_quad_rule.empty()) {
                auto fit = collisions.face_collisions.find(f);
                if (fit != collisions.face_collisions.end()) {
                    for (size_t qi = 0; qi < face_quad_rule.size(); qi++) {
                        if (qi >= fit->second.size())
                            continue;
                        const auto& qp = face_quad_rule[qi];
                        const Eigen::RowVector3d q_pos = qp.lambda[0]
                                * rest_positions.row(mesh.faces()(f, 0))
                            + qp.lambda[1]
                                * rest_positions.row(mesh.faces()(f, 1))
                            + qp.lambda[2]
                                * rest_positions.row(mesh.faces()(f, 2));
                        const auto& dict = *fit->second[qi];
                        for (int ci = 0; ci < dict.size(); ci++) {
                            auto pvids = get_primitive_vids(dict[ci]);
                            if (!pvids.empty()) {
                                active_pairs.push_back(
                                    { &dict[ci], true, q_pos,
                                      std::move(pvids) });
                            }
                        }
                    }
                }
            } else {
                // Vertex collisions when face quadrature is not used
                for (int lv = 0; lv < 3; lv++) {
                    const index_t v = mesh.faces()(f, lv);
                    auto vit = collisions.vertex_collisions.find(v);
                    if (vit == collisions.vertex_collisions.end()) {
                        continue;
                    }
                    const auto& dict = *vit->second;
                    for (int ci = 0; ci < dict.size(); ci++) {
                        auto pvids = get_primitive_vids(dict[ci]);
                        if (!pvids.empty()) {
                            active_pairs.push_back(
                                { &dict[ci], false, {}, std::move(pvids) });
                        }
                    }
                }
            }

            // Edge-edge collisions
            for (int le = 0; le < 3; le++) {
                const index_t edge_id = mesh.faces_to_edges()(f, le);
                const index_t ea = mesh.edges()(edge_id, 0);
                const index_t eb = mesh.edges()(edge_id, 1);

                for (index_t other_edge_id :
                     collisions.m_candidates.ee_set(edge_id)) {
                    const index_t ec = mesh.edges()(other_edge_id, 0);
                    const index_t ed = mesh.edges()(other_edge_id, 1);
                    if (ea == ec || ea == ed || eb == ec || eb == ed)
                        continue;

                    auto eit = collisions.edge_edge_collisions.find(
                        std::make_pair(edge_id, other_edge_id));
                    if (eit == collisions.edge_edge_collisions.end())
                        continue;

                    const auto& dict = *eit->second;
                    if (dict.ee_dtype() != EdgeEdgeDistanceType::EA_EB)
                        continue;

                    const double uv = closest_point_uv<double>(
                        rest_positions.row(ea), rest_positions.row(eb),
                        rest_positions.row(ec), rest_positions.row(ed),
                        dict.ee_dtype());
                    const Eigen::RowVector3d ee_qp =
                        uv * (rest_positions.row(eb) - rest_positions.row(ea))
                        + rest_positions.row(ea);

                    for (int ci = 0; ci < dict.size(); ci++) {
                        auto pvids = get_primitive_vids(dict[ci]);
                        if (!pvids.empty()) {
                            active_pairs.push_back(
                                { &dict[ci], true, ee_qp, std::move(pvids) });
                        }
                    }
                }
            }
        }

        std::vector<bool> completed(active_pairs.size(), false);
        bool has_active = true;
        while (has_active) {
            has_active = false;
            std::vector<bool> needs_reduction(nv, false);
            int num_remaining = 0;

            for (size_t i = 0; i < active_pairs.size(); i++) {
                if (completed[i])
                    continue;
                auto& p = active_pairs[i];
                const Eigen::VectorXd dofs = p.needs_extended
                    ? p.cc->dof(VertexMatrixView<3>(rest_positions, p.qp_pos))
                    : p.cc->dof(rest_positions);
                const double val = p.cc->weight * (*p.cc)(dofs, params, this);

                if (val != 0.0) {
                    for (const index_t vid : p.primitive_vids)
                        needs_reduction[vid] = true;
                    has_active = true;
                } else {
                    completed[i] = true;
                }
            }

            for (int v = 0; v < nv; v++) {
                if (needs_reduction[v])
                    m_values(v) *= zeta;
            }
        }

    } else if (mesh.dim() == 2) {
        struct ActivePair2D {
            const EspCollision* cc;
            bool needs_extended;
            Eigen::RowVector2d qp_pos;
            std::vector<index_t> primitive_vids;
        };

        std::vector<ActivePair2D> active_pairs;
        const GaussLobatto::Rule& rule =
            GaussLobatto::get_rule(params.quad_order);

        for (const auto& [ei, qp_dicts] : collisions.edge_collisions_2d) {
            const index_t e0 = mesh.edges()(ei, 0);
            const index_t e1 = mesh.edges()(ei, 1);
            for (size_t qi = 0; qi < qp_dicts.size(); qi++) {
                const auto& dict = *qp_dicts[qi];
                if (dict.size() == 0)
                    continue;
                const auto& qp = rule[qi];
                const Eigen::RowVector2d q_pos =
                    (1.0 - qp.xi) * rest_positions.row(e0)
                    + qp.xi * rest_positions.row(e1);
                for (int ci = 0; ci < dict.size(); ci++) {
                    auto pvids = get_primitive_vids(dict[ci]);
                    if (!pvids.empty())
                        active_pairs.push_back(
                            { &dict[ci], true, q_pos, std::move(pvids) });
                }
            }
        }

        std::vector<bool> completed(active_pairs.size(), false);
        bool has_active = true;
        while (has_active) {
            has_active = false;
            std::vector<bool> needs_reduction(nv, false);

            for (size_t i = 0; i < active_pairs.size(); i++) {
                if (completed[i])
                    continue;
                auto& p = active_pairs[i];
                const Eigen::VectorXd dofs = p.needs_extended
                    ? p.cc->dof(VertexMatrixView<2>(rest_positions, p.qp_pos))
                    : p.cc->dof(rest_positions);
                const double val = p.cc->weight * (*p.cc)(dofs, params, this);
                if (val != 0.0) {
                    for (const index_t vid : p.primitive_vids)
                        needs_reduction[vid] = true;
                    has_active = true;
                } else {
                    completed[i] = true;
                }
            }

            for (int v = 0; v < nv; v++) {
                if (needs_reduction[v])
                    m_values(v) *= zeta;
            }
        }
    }
}

double AdaptiveSupport::vertex(index_t vertex_id) const
{
    return m_values(vertex_id);
}

double AdaptiveSupport::edge(index_t edge_id, double t) const
{
    const int v0 = m_mesh->edges()(edge_id, 0);
    const int v1 = m_mesh->edges()(edge_id, 1);
    return (1.0 - t) * m_values(v0) + t * m_values(v1);
}

double AdaptiveSupport::face(index_t face_id, double u, double v) const
{
    const int f0 = m_mesh->faces()(face_id, 0);
    const int f1 = m_mesh->faces()(face_id, 1);
    const int f2 = m_mesh->faces()(face_id, 2);
    const double w = 1.0 - u - v;
    return w * m_values(f0) + u * m_values(f1) + v * m_values(f2);
}

} // namespace ipc
