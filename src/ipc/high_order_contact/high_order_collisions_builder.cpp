#include "high_order_collisions_builder.hpp"

#include "collisions/high_order_quadrature.hpp"

#include <ipc/distance/distance_type_exact.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/high_order_contact/quadrature_potential.hpp>

#include <tbb/enumerable_thread_specific.h>

#include <array>
#include <sstream>

namespace ipc {

using IntegrationType = HighOrderContactParameters::IntegrationType;

void HighOrderCollisionsBuilder<2>::build_edge_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const Candidates& candidates,
    const HighOrderContactParameters& params,
    size_t start,
    size_t end)
{
    const PointPotential pp(mesh, candidates, params);
    const GaussLobatto::Rule& rule = GaussLobatto::get_rule(params.quad_order);

    for (size_t edge_idx = start; edge_idx < end; ++edge_idx) {
        const index_t ei = static_cast<index_t>(edge_idx);

        if (candidates.ev_set(ei).empty() && candidates.ee_set(ei).empty())
            continue;

        if (params.integration_type == IntegrationType::NO_OBST
            && mesh.is_obstacle_edge(ei))
            continue;
        if (params.integration_type != IntegrationType::BRUTE_FORCE
            && mesh.is_obstacle_edge(ei)) {
            const auto& ev = candidates.ev_set(ei);
            const bool has_non_obstacle =
                std::any_of(ev.begin(), ev.end(), [&](index_t v) {
                    return !mesh.is_obstacle_vertex(v);
                });
            if (!has_non_obstacle)
                continue;
        }

        const double dhat = params.dhat;
        std::vector<std::unique_ptr<HighOrderCollisionDict<PointType::EDGE, 2>>>
            qp_dicts;
        qp_dicts.reserve(rule.size());
        bool has_any = false;

        for (const auto& qp : rule) {
            const std::array<double, 2> lambda = { { 1.0 - qp.xi, qp.xi } };
            size_t n = 0;
            auto dict = pp.build_collisions_at_edge_qp(V, ei, lambda, dhat, n);
            if (dict && dict->size() > 0)
                has_any = true;
            qp_dicts.push_back(std::move(dict));
        }

        if (has_any) {
            edge_collisions_2d.emplace_back(ei, std::move(qp_dicts));
        }
    }
}

void HighOrderCollisionsBuilder<2>::merge(
    tbb::enumerable_thread_specific<HighOrderCollisionsBuilder<2>>&
        local_storage,
    HighOrderCollisions& merged_collisions)
{
    size_t total_pairs = 0;

    // Move per-edge dicts into merged_collisions. No edge is processed by
    // more than one thread, so there are no duplicate edge keys.
    for (auto& builder : local_storage) {
        for (auto& [ei, dicts] : builder.edge_collisions_2d) {
            for (const auto& dict : dicts) {
                total_pairs += dict->size();
            }
            // Use insert with a value_type pair to force the move path.
            merged_collisions.edge_collisions_2d.insert(
                std::make_pair(ei, std::move(dicts)));
        }
    }

    logger().trace("2D edge QP collision pairs: {}.", total_pairs);
}

// ============================================================================

std::shared_ptr<HighOrderCollision>
HighOrderCollisionsBuilder<3>::reduce_point_triangle_collision(
    const FaceVertexCandidate& candidate,
    const HighOrderContactParameters& params,
    const CollisionMesh& mesh,
    const VertexMatrixView<3>& vertices,
    PointTriangleDistanceType dtype)
{
    const index_t vi = candidate.vertex_id;
    const index_t fi = candidate.face_id;

    const index_t t0 = mesh.faces()(fi, 0);
    const index_t t1 = mesh.faces()(fi, 1);
    const index_t t2 = mesh.faces()(fi, 2);

    const index_t e0 = mesh.faces_to_edges()(fi, 0);
    const index_t e1 = mesh.faces_to_edges()(fi, 1);
    const index_t e2 = mesh.faces_to_edges()(fi, 2);

    assert(vi != t0 && vi != t1 && vi != t2);

    if (dtype == PointTriangleDistanceType::AUTO) {
        dtype = point_triangle_distance_type_exact(
            vertices(vi), vertices(t0), vertices(t1), vertices(t2));
    }

    const double dist_sqr = point_triangle_distance(
        vertices(vi), vertices(t0), vertices(t1), vertices(t2), dtype);
    if (dist_sqr >= params.dhat * params.dhat) {
        return nullptr;
    }

    switch (dtype) {
    case PointTriangleDistanceType::P_T0:
        return std::make_shared<HighOrderCollisionTemplate<Vertex3, Vertex3>>(
            t0, vi, mesh);

    case PointTriangleDistanceType::P_T1:
        return std::make_shared<HighOrderCollisionTemplate<Vertex3, Vertex3>>(
            t1, vi, mesh);

    case PointTriangleDistanceType::P_T2:
        return std::make_shared<HighOrderCollisionTemplate<Vertex3, Vertex3>>(
            t2, vi, mesh);

    case PointTriangleDistanceType::P_E0:
        return std::make_shared<HighOrderCollisionTemplate<Edge3P1, Vertex3>>(
            e0, vi, mesh);

    case PointTriangleDistanceType::P_E1:
        return std::make_shared<HighOrderCollisionTemplate<Edge3P1, Vertex3>>(
            e1, vi, mesh);

    case PointTriangleDistanceType::P_E2:
        return std::make_shared<HighOrderCollisionTemplate<Edge3P1, Vertex3>>(
            e2, vi, mesh);

    case PointTriangleDistanceType::P_T:
        return std::make_shared<HighOrderCollisionTemplate<Face3P1, Vertex3>>(
            fi, vi, mesh);

    case PointTriangleDistanceType::AUTO:
    default:
        assert(false);
        return std::make_shared<HighOrderCollisionTemplate<Face3P1, Vertex3>>(
            fi, vi, mesh);
    }
}

std::shared_ptr<HighOrderCollision>
HighOrderCollisionsBuilder<3>::reduce_point_edge_collision(
    const EdgeVertexCandidate& candidate,
    const HighOrderContactParameters& params,
    const CollisionMesh& mesh,
    const VertexMatrixView<3>& vertices,
    PointEdgeDistanceType dtype)
{
    const index_t vi = candidate.vertex_id;
    const index_t ei = candidate.edge_id;

    const index_t t0 = mesh.edges()(ei, 0);
    const index_t t1 = mesh.edges()(ei, 1);

    if (dtype == PointEdgeDistanceType::AUTO) {
        dtype = point_edge_distance_type_exact(
            vertices(vi), vertices(t0), vertices(t1));
    }

    const double dist_sqr =
        point_edge_distance(vertices(vi), vertices(t0), vertices(t1), dtype);
    if (dist_sqr >= params.dhat * params.dhat) {
        return nullptr;
    }

    switch (dtype) {
    case PointEdgeDistanceType::P_E0:
        return std::make_shared<HighOrderCollisionTemplate<Vertex3, Vertex3>>(
            t0, vi, mesh);
    case PointEdgeDistanceType::P_E1:
        return std::make_shared<HighOrderCollisionTemplate<Vertex3, Vertex3>>(
            t1, vi, mesh);
    case PointEdgeDistanceType::P_E:
        return std::make_shared<HighOrderCollisionTemplate<Edge3P1, Vertex3>>(
            ei, vi, mesh);
    default:
        assert(false);
        return std::make_shared<HighOrderCollisionTemplate<Edge3P1, Vertex3>>(
            ei, vi, mesh);
    }
}

// ============================================================================
// QuadratureCollisionsBuilder

QuadratureCollisionsBuilder::QuadratureCollisionsBuilder(
    const CollisionMesh& mesh,
    const Candidates& candidates,
    const HighOrderContactParameters& params)
    : point_potential(
          std::make_shared<PointPotential>(mesh, candidates, params))
{
}

QuadratureCollisionsBuilder::~QuadratureCollisionsBuilder() = default;

QuadratureCollisionsBuilder::QuadratureCollisionsBuilder(
    const QuadratureCollisionsBuilder& other)
{
    point_potential = other.point_potential;
    vertex_collisions.clear();
    for (const auto& cc : other.vertex_collisions) {
        vertex_collisions.push_back(
            std::make_unique<HighOrderCollisionDict<PointType::VERTEX>>(*cc));
    }
    edge_edge_collisions.clear();
    for (const auto& cc : other.edge_edge_collisions) {
        edge_edge_collisions.push_back(
            std::make_unique<HighOrderCollisionDict<PointType::EDGE>>(*cc));
    }
    face_collisions.clear();
    for (const auto& [fi, dicts] : other.face_collisions) {
        std::vector<std::unique_ptr<HighOrderCollisionDict<PointType::FACE>>>
            copied;
        for (const auto& d : dicts) {
            copied.push_back(
                std::make_unique<HighOrderCollisionDict<PointType::FACE>>(*d));
        }
        face_collisions.push_back({ fi, std::move(copied) });
    }
}
QuadratureCollisionsBuilder&
QuadratureCollisionsBuilder::operator=(const QuadratureCollisionsBuilder& other)
{
    point_potential = other.point_potential;
    vertex_collisions.clear();
    for (const auto& cc : other.vertex_collisions) {
        vertex_collisions.push_back(
            std::make_unique<HighOrderCollisionDict<PointType::VERTEX>>(*cc));
    }
    edge_edge_collisions.clear();
    for (const auto& cc : other.edge_edge_collisions) {
        edge_edge_collisions.push_back(
            std::make_unique<HighOrderCollisionDict<PointType::EDGE>>(*cc));
    }
    face_collisions.clear();
    for (const auto& [fi, dicts] : other.face_collisions) {
        std::vector<std::unique_ptr<HighOrderCollisionDict<PointType::FACE>>>
            copied;
        for (const auto& d : dicts) {
            copied.push_back(
                std::make_unique<HighOrderCollisionDict<PointType::FACE>>(*d));
        }
        face_collisions.push_back({ fi, std::move(copied) });
    }
    return *this;
}

void QuadratureCollisionsBuilder::build_vertex_collisions(
    const Eigen::MatrixXd& vertices,
    const std::vector<index_t>& vertex_indices,
    const size_t start_i,
    const size_t end_i)
{
    const CollisionMesh& mesh = point_potential->mesh;
    const HighOrderContactParameters& params = point_potential->params;
    for (size_t i = start_i; i < end_i; i++) {
        const index_t vi = vertex_indices[i];
        if (params.integration_type == IntegrationType::NO_OBST
            && mesh.is_obstacle_vertex(vi))
            continue;
        if (params.integration_type != IntegrationType::BRUTE_FORCE
            && mesh.is_obstacle_vertex(vi)) {
            const auto v_set = point_potential->candidates.vv_set(vi);
            const auto e_set = point_potential->candidates.ve_set(vi);
            const auto f_set = point_potential->candidates.vf_set(vi);
            const bool has_non_obstacle =
                std::any_of(
                    v_set.begin(), v_set.end(),
                    [&](index_t v) { return !mesh.is_obstacle_vertex(v); })
                || std::any_of(
                    e_set.begin(), e_set.end(),
                    [&](index_t e) { return !mesh.is_obstacle_edge(e); })
                || std::any_of(f_set.begin(), f_set.end(), [&](index_t f) {
                       return !mesh.is_obstacle_face(f);
                   });
            if (!has_non_obstacle)
                continue;
        }
        size_t n = 0;
        auto dict =
            point_potential->build_collisions_at_vertex(vertices, vi, n);
        if (dict && dict->size() > 0) {
            vertex_collisions.push_back(std::move(dict));
        }
        num_collision_pairs += n;
    }
}

void QuadratureCollisionsBuilder::build_face_collisions(
    const Eigen::MatrixXd& vertices,
    const std::vector<index_t>& face_indices,
    const size_t start_i,
    const size_t end_i)
{
    const CollisionMesh& mesh = point_potential->mesh;
    const auto& face_quad_rule = point_potential->params.get_quad_rule();
    if (face_quad_rule.empty())
        return;

    const HighOrderContactParameters& params = point_potential->params;
    for (size_t i = start_i; i < end_i; i++) {
        const index_t fi = face_indices[i];
        if (params.integration_type == IntegrationType::NO_OBST
            && mesh.is_obstacle_face(fi))
            continue;
        if (params.integration_type != IntegrationType::BRUTE_FORCE
            && mesh.is_obstacle_face(fi)) {
            const auto v_set = point_potential->candidates.fv_set(fi);
            const auto e_set = point_potential->candidates.fe_set(fi);
            const auto f_set = point_potential->candidates.ff_set(fi);
            const bool has_non_obstacle =
                std::any_of(
                    v_set.begin(), v_set.end(),
                    [&](index_t v) { return !mesh.is_obstacle_vertex(v); })
                || std::any_of(
                    e_set.begin(), e_set.end(),
                    [&](index_t e) { return !mesh.is_obstacle_edge(e); })
                || std::any_of(f_set.begin(), f_set.end(), [&](index_t f) {
                       return !mesh.is_obstacle_face(f);
                   });
            if (!has_non_obstacle)
                continue;
        }

        std::vector<std::unique_ptr<HighOrderCollisionDict<PointType::FACE>>>
            per_qp_dicts;
        per_qp_dicts.reserve(face_quad_rule.size());
        bool any_nonempty = false;
        for (const auto& qp : face_quad_rule) {
            size_t n = 0;
            auto dict =
                point_potential->build_collisions_at_face_interior_point(
                    vertices, fi, qp.lambda, n);
            if (dict && dict->size() > 0) {
                any_nonempty = true;
            }
            // Keep the qi indexing aligned with face_quad_rule, even if the
            // dict is empty for this quadrature point.
            per_qp_dicts.push_back(std::move(dict));
            num_collision_pairs += n;
        }
        if (any_nonempty) {
            face_collisions.push_back({ fi, std::move(per_qp_dicts) });
        }
    }
}

void QuadratureCollisionsBuilder::build_edge_edge_collisions(
    const Eigen::MatrixXd& vertices,
    const std::vector<EdgeEdgeCandidate>& ee_candidates,
    const size_t start_i,
    const size_t end_i)
{
    const HighOrderContactParameters& params = point_potential->params;
    const CollisionMesh& mesh = point_potential->mesh;

    // Returns true if edge e (which is an obstacle) has at least one
    // non-obstacle candidate. Used in NORMAL mode to skip placing a QP on an
    // obstacle edge with only obstacle candidates.
    auto obstacle_edge_has_non_obstacle_candidates = [&](index_t e) -> bool {
        const auto v_set = point_potential->candidates.ev_set(e);
        const auto e_set = point_potential->candidates.ee_set(e);
        const auto f_set = point_potential->candidates.ef_set(e);
        return std::any_of(
                   v_set.begin(), v_set.end(),
                   [&](index_t v) { return !mesh.is_obstacle_vertex(v); })
            || std::any_of(
                   e_set.begin(), e_set.end(),
                   [&](index_t e2) { return !mesh.is_obstacle_edge(e2); })
            || std::any_of(f_set.begin(), f_set.end(), [&](index_t f) {
                   return !mesh.is_obstacle_face(f);
               });
    };

    for (size_t i = start_i; i < end_i; i++) {
        const auto& candidate = ee_candidates[i];
        const index_t ei = candidate.edge0_id;
        const index_t ej = candidate.edge1_id;

        const index_t ea = mesh.edges()(ei, 0);
        const index_t eb = mesh.edges()(ei, 1);
        const index_t ec = mesh.edges()(ej, 0);
        const index_t ed = mesh.edges()(ej, 1);

        if (ea == ec || ea == ed || eb == ec || eb == ed) {
            continue;
        }

        if (params.integration_type != IntegrationType::BRUTE_FORCE
            && mesh.is_obstacle_edge(ei) && mesh.is_obstacle_edge(ej)) {
            continue;
        }

        if (is_parallel_edge_edge(
                vertices.row(ea), vertices.row(eb), vertices.row(ec),
                vertices.row(ed))) {
            continue;
        }

        const auto dtype = edge_edge_distance_type_exact(
            vertices.row(ea), vertices.row(eb), vertices.row(ec),
            vertices.row(ed));

        const double dist_sq = edge_edge_distance(
            vertices.row(ea), vertices.row(eb), vertices.row(ec),
            vertices.row(ed), dtype);

        if (dist_sq >= params.dbar * params.dbar) {
            continue;
        }

        // HighOrderContactPotential only ever evaluates dicts whose stored
        // dtype is EA_EB (see the `if (dtype != EA_EB) continue;` guards in
        // high_order_contact_potential.cpp). All other edge-edge distance
        // types are captured through vertex_collisions at the relevant
        // endpoint, so building EA_EB0/EA_EB1/EA0_EB/EA1_EB dicts here is
        // dead work.
        if (dtype != EdgeEdgeDistanceType::EA_EB) {
            continue;
        }

        const bool ei_is_obs = mesh.is_obstacle_edge(ei);
        const bool ej_is_obs = mesh.is_obstacle_edge(ej);

        if ((params.integration_type != IntegrationType::NO_OBST || !ei_is_obs)
            && (!ei_is_obs
                || params.integration_type == IntegrationType::BRUTE_FORCE
                || obstacle_edge_has_non_obstacle_candidates(ei))) {
            size_t n = 0;
            auto dict =
                point_potential->build_collisions_at_edge_edge_closest_point(
                    vertices, ei, ej, dtype, n);
            if (dict && dict->size() > 0) {
                edge_edge_collisions.push_back(std::move(dict));
            }
            num_collision_pairs += n;
        }

        if ((params.integration_type != IntegrationType::NO_OBST || !ej_is_obs)
            && (!ej_is_obs
                || params.integration_type == IntegrationType::BRUTE_FORCE
                || obstacle_edge_has_non_obstacle_candidates(ej))) {
            size_t n = 0;
            auto dict =
                point_potential->build_collisions_at_edge_edge_closest_point(
                    vertices, ej, ei, dtype, n);
            if (dict && dict->size() > 0) {
                edge_edge_collisions.push_back(std::move(dict));
            }
            num_collision_pairs += n;
        }
    }
}

void QuadratureCollisionsBuilder::merge(
    tbb::enumerable_thread_specific<QuadratureCollisionsBuilder>& local_storage,
    HighOrderCollisions& merged_collisions)
{
    // Reserve space
    size_t total_v = 0, total_ee = 0, total_f = 0;
    for (const auto& storage : local_storage) {
        total_v += storage.vertex_collisions.size();
        total_ee += storage.edge_edge_collisions.size();
        total_f += storage.face_collisions.size();
    }
    merged_collisions.vertex_collisions.reserve(total_v);
    merged_collisions.edge_edge_collisions.reserve(total_ee);
    merged_collisions.face_collisions.reserve(total_f);

    for (auto& storage : local_storage) {
        for (auto& cc : storage.vertex_collisions) {
            merged_collisions.vertex_collisions.insert(
                std::make_pair<
                    index_t,
                    std::unique_ptr<HighOrderCollisionDict<PointType::VERTEX>>>(
                    cc->primitive_id(), std::move(cc)));
        }
        for (auto& cc : storage.edge_edge_collisions) {
            const auto id = cc->primitive_ids();
            merged_collisions.edge_edge_collisions.insert(
                std::make_pair(std::make_pair(id[0], id[1]), std::move(cc)));
        }
        for (auto& [fi, dicts] : storage.face_collisions) {
            merged_collisions.face_collisions.emplace(fi, std::move(dicts));
        }
        merged_collisions.num_quadrature_collision_pairs +=
            storage.num_collision_pairs;
    }
}

} // namespace ipc
