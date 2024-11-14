#include "normal_collisions_builder.hpp"

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_triangle.hpp>

namespace ipc {

NormalCollisionsBuilder::NormalCollisionsBuilder(
    const bool _use_area_weighting, const bool _enable_shape_derivatives)
    : use_area_weighting(_use_area_weighting)
    , enable_shape_derivatives(_enable_shape_derivatives)
{
}

// ============================================================================

void NormalCollisionsBuilder::add_vertex_vertex_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<VertexVertexCandidate>& candidates,
    const std::function<bool(double)>& is_active,
    const size_t start_i,
    const size_t end_i)
{
    for (size_t i = start_i; i < end_i; i++) {
        const auto& [vi, vj] = candidates[i];

        const double distance =
            point_point_distance(vertices.row(vi), vertices.row(vj));
        point_point_distance(vertices.row(vi), vertices.row(vj));
        if (!is_active(distance)) {
            continue;
        }

        // ÷ 2 to handle double counting. Sum vertex areas because duplicate
        // vertex-vertex candidates were removed.
        const double weight = use_area_weighting
            ? (0.5 * (mesh.vertex_area(vi) + mesh.vertex_area(vj)))
            : 1;

        Eigen::SparseVector<double> weight_gradient;
        if (enable_shape_derivatives) {
            weight_gradient = use_area_weighting
                ? (0.5
                   * (mesh.vertex_area_gradient(vi)
                      + mesh.vertex_area_gradient(vj)))
                : Eigen::SparseVector<double>(vertices.size());
        }

        VertexVertexNormalCollision vv(vi, vj, weight, weight_gradient);
        vv_to_id.emplace(vv, vv_collisions.size());
        vv_collisions.push_back(vv);
    }
}

void NormalCollisionsBuilder::add_edge_vertex_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<EdgeVertexCandidate>& candidates,
    const std::function<bool(double)>& is_active,
    const size_t start_i,
    const size_t end_i)
{
    for (size_t i = start_i; i < end_i; i++) {
        const auto& [ei, vi] = candidates[i];
        const auto [v, e0, e1, _] =
            candidates[i].vertices(vertices, mesh.edges(), mesh.faces());
        const PointEdgeDistanceType dtype = point_edge_distance_type(v, e0, e1);
        const double distance_sqr = point_edge_distance(v, e0, e1, dtype);

        if (!is_active(distance_sqr))
            continue;

        // ÷ 2 to handle double counting for correct integration
        const double weight =
            use_area_weighting ? (0.5 * mesh.vertex_area(vi)) : 1;

        Eigen::SparseVector<double> weight_gradient;
        if (enable_shape_derivatives) {
            weight_gradient = use_area_weighting
                ? (0.5 * mesh.vertex_area_gradient(vi))
                : Eigen::SparseVector<double>(vertices.size());
        }

        add_edge_vertex_collision(
            mesh, candidates[i], dtype, weight, weight_gradient);
    }
}

void NormalCollisionsBuilder::add_edge_vertex_collision(
    const CollisionMesh& mesh,
    const EdgeVertexCandidate& candidate,
    const PointEdgeDistanceType dtype,
    const double weight,
    const Eigen::SparseVector<double>& weight_gradient)
{
    const auto& [ei, vi] = candidate;

    switch (dtype) {
    case PointEdgeDistanceType::P_E0:
        add_vertex_vertex_collision(
            vi, mesh.edges()(ei, 0), weight, weight_gradient);
        break;

    case PointEdgeDistanceType::P_E1:
        add_vertex_vertex_collision(
            vi, mesh.edges()(ei, 1), weight, weight_gradient);
        break;

    case PointEdgeDistanceType::P_E:
        add_edge_vertex_collision(ei, vi, weight, weight_gradient);
        break;

    case PointEdgeDistanceType::AUTO:
    default:
        assert(false);
        break;
    }
}

void NormalCollisionsBuilder::add_edge_edge_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<EdgeEdgeCandidate>& candidates,
    const std::function<bool(double)>& is_active,
    const size_t start_i,
    const size_t end_i)
{
    for (size_t i = start_i; i < end_i; i++) {
        const auto& [eai, ebi] = candidates[i];

        const auto [ea0i, ea1i, eb0i, eb1i] =
            candidates[i].vertex_ids(mesh.edges(), mesh.faces());

        const auto [ea0, ea1, eb0, eb1] =
            candidates[i].vertices(vertices, mesh.edges(), mesh.faces());

        const EdgeEdgeDistanceType actual_dtype =
            edge_edge_distance_type(ea0, ea1, eb0, eb1);

        const double distance_sqr =
            edge_edge_distance(ea0, ea1, eb0, eb1, actual_dtype);

        if (!is_active(distance_sqr))
            continue;

        const double eps_x = edge_edge_mollifier_threshold(
            mesh.rest_positions().row(ea0i), mesh.rest_positions().row(ea1i),
            mesh.rest_positions().row(eb0i), mesh.rest_positions().row(eb1i));

        const double ee_cross_norm_sqr =
            edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1);

        // NOTE: This may not actually be the distance type, but all EE
        // pairs requiring mollification must be mollified later.
        const EdgeEdgeDistanceType dtype = ee_cross_norm_sqr < eps_x
            ? EdgeEdgeDistanceType::EA_EB
            : actual_dtype;

        // ÷ 4 to handle double counting and PT + EE for correct integration.
        // Sum edge areas because duplicate edge candidates were removed.
        const double weight = use_area_weighting
            ? (0.25 * (mesh.edge_area(eai) + mesh.edge_area(ebi)))
            : 1;

        Eigen::SparseVector<double> weight_gradient;
        if (enable_shape_derivatives) {
            weight_gradient = use_area_weighting
                ? (0.25
                   * (mesh.edge_area_gradient(eai)
                      + mesh.edge_area_gradient(ebi)))
                : Eigen::SparseVector<double>(vertices.size());
        }

        switch (dtype) {
        case EdgeEdgeDistanceType::EA0_EB0:
            add_vertex_vertex_collision(ea0i, eb0i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA0_EB1:
            add_vertex_vertex_collision(ea0i, eb1i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA1_EB0:
            add_vertex_vertex_collision(ea1i, eb0i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA1_EB1:
            add_vertex_vertex_collision(ea1i, eb1i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA_EB0:
            add_edge_vertex_collision(eai, eb0i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA_EB1:
            add_edge_vertex_collision(eai, eb1i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA0_EB:
            add_edge_vertex_collision(ebi, ea0i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA1_EB:
            add_edge_vertex_collision(ebi, ea1i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA_EB:
            ee_collisions.emplace_back(
                eai, ebi, eps_x, weight, weight_gradient, actual_dtype);
            ee_to_id.emplace(ee_collisions.back(), ee_collisions.size() - 1);
            break;

        case EdgeEdgeDistanceType::AUTO:
        default:
            assert(false);
            break;
        }
    }
}

void NormalCollisionsBuilder::add_face_vertex_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<FaceVertexCandidate>& candidates,
    const std::function<bool(double)>& is_active,
    const size_t start_i,
    const size_t end_i)
{
    for (size_t i = start_i; i < end_i; i++) {
        const auto& [fi, vi] = candidates[i];
        const long f0i = mesh.faces()(fi, 0), f1i = mesh.faces()(fi, 1),
                   f2i = mesh.faces()(fi, 2);

        const auto [v, f0, f1, f2] =
            candidates[i].vertices(vertices, mesh.edges(), mesh.faces());

        // Compute distance type
        const PointTriangleDistanceType dtype =
            point_triangle_distance_type(v, f0, f1, f2);
        const double distance_sqr =
            point_triangle_distance(v, f0, f1, f2, dtype);

        if (!is_active(distance_sqr))
            continue;

        // ÷ 4 to handle double counting and PT + EE for correct integration
        const double weight =
            use_area_weighting ? (0.25 * mesh.vertex_area(vi)) : 1;

        Eigen::SparseVector<double> weight_gradient;
        if (enable_shape_derivatives) {
            weight_gradient = use_area_weighting
                ? (0.25 * mesh.vertex_area_gradient(vi))
                : Eigen::SparseVector<double>(vertices.size());
        }

        switch (dtype) {
        case PointTriangleDistanceType::P_T0:
            add_vertex_vertex_collision(vi, f0i, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::P_T1:
            add_vertex_vertex_collision(vi, f1i, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::P_T2:
            add_vertex_vertex_collision(vi, f2i, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::P_E0:
            add_edge_vertex_collision(
                mesh.faces_to_edges()(fi, 0), vi, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::P_E1:
            add_edge_vertex_collision(
                mesh.faces_to_edges()(fi, 1), vi, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::P_E2:
            add_edge_vertex_collision(
                mesh.faces_to_edges()(fi, 2), vi, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::P_T:
            fv_collisions.emplace_back(fi, vi, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::AUTO:
        default:
            assert(false);
            break;
        }
    }
}

// ============================================================================

void NormalCollisionsBuilder::add_edge_vertex_negative_vertex_vertex_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<VertexVertexCandidate>& candidates,
    const size_t start_i,
    const size_t end_i)
{
    const auto add_weight = [&](const size_t vi, const size_t vj,
                                double& weight,
                                Eigen::SparseVector<double>& weight_gradient) {
        const auto& incident_vertices = mesh.vertex_vertex_adjacencies()[vj];
        const int incident_edge_amt = incident_vertices.size()
            - int(incident_vertices.find(vi) != incident_vertices.end());

        if (incident_edge_amt > 1) {
            // ÷ 2 to handle double counting for correct integration
            weight += (1 - incident_edge_amt)
                * (use_area_weighting ? (0.5 * mesh.vertex_area(vi)) : 1);

            if (enable_shape_derivatives && use_area_weighting) {
                weight_gradient += 0.5 * (1 - incident_edge_amt)
                    * mesh.vertex_area_gradient(vi);
            }
        }
    };

    for (size_t i = start_i; i < end_i; i++) {
        const auto& [vi, vj] = candidates[i];
        assert(vi != vj);

        double weight = 0;
        Eigen::SparseVector<double> weight_gradient;
        if (enable_shape_derivatives) {
            weight_gradient = Eigen::SparseVector<double>(vertices.size());
        }

        add_weight(vi, vj, weight, weight_gradient);
        add_weight(vj, vi, weight, weight_gradient);

        if (weight != 0) {
            add_vertex_vertex_collision(vi, vj, weight, weight_gradient);
        }
    }
}

void NormalCollisionsBuilder::add_face_vertex_positive_vertex_vertex_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<VertexVertexCandidate>& candidates,
    const size_t start_i,
    const size_t end_i)
{
    const auto add_weight = [&](const size_t vi, const size_t vj,
                                double& weight,
                                Eigen::SparseVector<double>& weight_gradient) {
        const auto& incident_vertices = mesh.vertex_vertex_adjacencies()[vj];
        if (mesh.is_vertex_on_boundary(vj)
            || incident_vertices.find(vi) != incident_vertices.end()) {
            return; // Skip boundary vertices and incident vertices
        }

        // ÷ 4 to handle double counting and PT + EE for correct integration.
        weight += use_area_weighting ? (0.25 * mesh.vertex_area(vi)) : 1;

        if (enable_shape_derivatives && use_area_weighting) {
            weight_gradient += 0.25 * mesh.vertex_area_gradient(vi);
        }
    };

    for (size_t i = start_i; i < end_i; i++) {
        const auto& [vi, vj] = candidates[i];
        assert(vi != vj);

        double weight = 0;
        Eigen::SparseVector<double> weight_gradient;
        if (enable_shape_derivatives) {
            weight_gradient = Eigen::SparseVector<double>(vertices.size());
        }

        add_weight(vi, vj, weight, weight_gradient);
        add_weight(vj, vi, weight, weight_gradient);

        if (weight != 0) {
            add_vertex_vertex_collision(vi, vj, weight, weight_gradient);
        }
    }
}

void NormalCollisionsBuilder::add_face_vertex_negative_edge_vertex_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<EdgeVertexCandidate>& candidates,
    const size_t start_i,
    const size_t end_i)
{
    for (size_t i = start_i; i < end_i; i++) {
        const auto& [ei, vi] = candidates[i];
        assert(vi != mesh.edges()(ei, 0) && vi != mesh.edges()(ei, 1));

        const auto& incident_vertices = mesh.edge_vertex_adjacencies()[ei];
        const int incident_triangle_amt = incident_vertices.size()
            - int(incident_vertices.find(vi) != incident_vertices.end());

        if (incident_triangle_amt > 1) {
            // ÷ 4 to handle double counting and PT + EE for correct integration
            const double weight = (1 - incident_triangle_amt)
                * (use_area_weighting ? (0.25 * mesh.vertex_area(vi)) : 1);

            Eigen::SparseVector<double> weight_gradient;
            if (enable_shape_derivatives) {
                weight_gradient = use_area_weighting
                    ? (0.25 * (1 - incident_triangle_amt)
                       * mesh.vertex_area_gradient(vi))
                    : Eigen::SparseVector<double>(vertices.size());
            }

            add_edge_vertex_collision(
                mesh, candidates[i],
                point_edge_distance_type(
                    vertices.row(vi), vertices.row(mesh.edges()(ei, 0)),
                    vertices.row(mesh.edges()(ei, 1))),
                weight, weight_gradient);
        }
    }
}

void NormalCollisionsBuilder::add_edge_edge_negative_edge_vertex_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<EdgeVertexCandidate>& candidates,
    const size_t start_i,
    const size_t end_i)
{
    // Notation: (ea, p) ∈ C, ea = (ea0, ea1) ∈ E, p ∈ eb = (p, q) ∈ E

    for (size_t i = start_i; i < end_i; i++) {
        const auto& [ea, p] = candidates[i];
        const int ea0 = mesh.edges()(ea, 0), ea1 = mesh.edges()(ea, 1);
        assert(p != ea0 && p != ea1);

        // ÷ 4 to handle double counting and PT + EE for correct integration
        const double weight =
            use_area_weighting ? (-0.25 * mesh.edge_area(ea)) : -1;
        Eigen::SparseVector<double> weight_gradient;
        if (enable_shape_derivatives) {
            weight_gradient = use_area_weighting
                ? (-0.25 * mesh.edge_area_gradient(ea))
                : Eigen::SparseVector<double>(vertices.size());
        }

        const PointEdgeDistanceType dtype = point_edge_distance_type(
            vertices.row(p), vertices.row(ea0), vertices.row(ea1));

        int nonmollified_incident_edge_amt = 0;

        const auto& incident_edges = mesh.vertex_edge_adjacencies()[p];
        for (const int eb : incident_edges) {
            const int eb0 = mesh.edges()(eb, 0), eb1 = mesh.edges()(eb, 1);
            const int q = mesh.edges()(eb, int(p == eb0));
            assert(p != q);
            if (q == ea0 || q == ea1) {
                continue;
            }

            const double eps_x = edge_edge_mollifier_threshold(
                mesh.rest_positions().row(ea0), mesh.rest_positions().row(ea1),
                mesh.rest_positions().row(eb0), mesh.rest_positions().row(eb1));

            const double ee_cross_norm_sqr = edge_edge_cross_squarednorm(
                vertices.row(ea0), vertices.row(ea1), vertices.row(eb0),
                vertices.row(eb1));

            if (ee_cross_norm_sqr >= eps_x) {
                nonmollified_incident_edge_amt++;
                continue;
            }

            // Add mollified EE collision with specified distance type
            // Convert the PE distance type to an EE distance type
            EdgeEdgeDistanceType ee_dtype = EdgeEdgeDistanceType::AUTO;
            switch (dtype) {
            case PointEdgeDistanceType::P_E0:
                ee_dtype = p == eb0 ? EdgeEdgeDistanceType::EA0_EB0
                                    : EdgeEdgeDistanceType::EA0_EB1;
                break;
            case PointEdgeDistanceType::P_E1:
                ee_dtype = p == eb0 ? EdgeEdgeDistanceType::EA1_EB0
                                    : EdgeEdgeDistanceType::EA1_EB1;
                break;
            case PointEdgeDistanceType::P_E:
                ee_dtype = p == eb0 ? EdgeEdgeDistanceType::EA_EB0
                                    : EdgeEdgeDistanceType::EA_EB1;
                break;
            default:
                assert(false);
                break;
            }

            add_edge_edge_collision(
                ea, eb, eps_x, weight, weight_gradient, ee_dtype);
        }

        if (nonmollified_incident_edge_amt == 1) {
            continue; // no collision to add because (ρ(x) - 1) = 0
        }
        // if nonmollified_incident_edge_amt == 0, then we need to explicitly
        // add a positive collision.
        add_edge_vertex_collision(
            mesh, candidates[i], dtype,
            (nonmollified_incident_edge_amt - 1) * weight,
            (nonmollified_incident_edge_amt - 1) * weight_gradient);
    }
}

// ============================================================================

void NormalCollisionsBuilder::add_vertex_vertex_collision(
    const VertexVertexNormalCollision& vv_collision,
    unordered_map<VertexVertexNormalCollision, long>& vv_to_id,
    std::vector<VertexVertexNormalCollision>& vv_collisions)
{
    auto found_item = vv_to_id.find(vv_collision);
    if (found_item != vv_to_id.end()) {
        // collision already exists, so increase weight
        vv_collisions[found_item->second].weight += vv_collision.weight;
        vv_collisions[found_item->second].weight_gradient +=
            vv_collision.weight_gradient;
    } else {
        // New collision, so add it to the end of vv_collisions
        vv_to_id.emplace(vv_collision, vv_collisions.size());
        vv_collisions.push_back(vv_collision);
    }
}

void NormalCollisionsBuilder::add_edge_vertex_collision(
    const EdgeVertexNormalCollision& ev_collision,
    unordered_map<EdgeVertexNormalCollision, long>& ev_to_id,
    std::vector<EdgeVertexNormalCollision>& ev_collisions)
{
    auto found_item = ev_to_id.find(ev_collision);
    if (found_item != ev_to_id.end()) {
        // collision already exists, so increase weight
        ev_collisions[found_item->second].weight += ev_collision.weight;
        ev_collisions[found_item->second].weight_gradient +=
            ev_collision.weight_gradient;
    } else {
        // New collision, so add it to the end of ev_collisions
        ev_to_id.emplace(ev_collision, ev_collisions.size());
        ev_collisions.push_back(ev_collision);
    }
}

void NormalCollisionsBuilder::add_edge_edge_collision(
    const EdgeEdgeNormalCollision& ee_collision,
    unordered_map<EdgeEdgeNormalCollision, long>& ee_to_id,
    std::vector<EdgeEdgeNormalCollision>& ee_collisions)
{
    auto found_item = ee_to_id.find(ee_collision);
    if (found_item != ee_to_id.end()) {
        // collision already exists, so increase weight
        assert(ee_collision == ee_collisions[found_item->second]);
        ee_collisions[found_item->second].weight += ee_collision.weight;
        ee_collisions[found_item->second].weight_gradient +=
            ee_collision.weight_gradient;
    } else {
        // New collision, so add it to the end of ee_collisions
        ee_to_id.emplace(ee_collision, ee_collisions.size());
        ee_collisions.push_back(ee_collision);
    }
}

// ============================================================================

void NormalCollisionsBuilder::merge(
    const tbb::enumerable_thread_specific<NormalCollisionsBuilder>&
        local_storage,
    NormalCollisions& merged_collisions)
{
    unordered_map<VertexVertexNormalCollision, long> vv_to_id;
    unordered_map<EdgeVertexNormalCollision, long> ev_to_id;
    unordered_map<EdgeEdgeNormalCollision, long> ee_to_id;
    auto& vv_collisions = merged_collisions.vv_collisions;
    auto& ev_collisions = merged_collisions.ev_collisions;
    auto& ee_collisions = merged_collisions.ee_collisions;
    auto& fv_collisions = merged_collisions.fv_collisions;

    // size up the hash items
    size_t n_vv = 0, n_ev = 0, n_ee = 0, n_fv = 0;
    for (const auto& storage : local_storage) {
        // This is an conservative estimate
        n_vv += storage.vv_collisions.size();
        n_ev += storage.ev_collisions.size();
        n_ee += storage.ee_collisions.size();
        n_fv += storage.fv_collisions.size();
    }
    vv_collisions.reserve(n_vv);
    ev_collisions.reserve(n_ev);
    ee_collisions.reserve(n_ee);
    fv_collisions.reserve(n_fv);

    // merge
    for (const auto& builder : local_storage) {
        if (vv_collisions.empty()) {
            vv_to_id = builder.vv_to_id;
            vv_collisions = builder.vv_collisions;
        } else {
            for (const auto& vv : builder.vv_collisions) {
                add_vertex_vertex_collision(vv, vv_to_id, vv_collisions);
            }
        }

        if (ev_collisions.empty()) {
            ev_to_id = builder.ev_to_id;
            ev_collisions = builder.ev_collisions;
        } else {
            for (const auto& ev : builder.ev_collisions) {
                add_edge_vertex_collision(ev, ev_to_id, ev_collisions);
            }
        }

        if (ee_collisions.empty()) {
            ee_to_id = builder.ee_to_id;
            ee_collisions = builder.ee_collisions;
        } else {
            for (const auto& ee : builder.ee_collisions) {
                add_edge_edge_collision(ee, ee_to_id, ee_collisions);
            }
        }

        fv_collisions.insert(
            fv_collisions.end(), builder.fv_collisions.begin(),
            builder.fv_collisions.end());
    }

    // If positive and negative vertex-vertex collisions cancel out, remove
    // them. This can happen when edge-vertex collisions reduce to
    // vertex-vertex collisions. This will avoid unnecessary computation.
    vv_collisions.erase(
        std::remove_if(
            vv_collisions.begin(), vv_collisions.end(),
            [&](const VertexVertexNormalCollision& vv) {
                return vv.weight == 0;
            }),
        vv_collisions.end());
    // Same for edge-vertex collisions.
    ev_collisions.erase(
        std::remove_if(
            ev_collisions.begin(), ev_collisions.end(),
            [&](const EdgeVertexNormalCollision& ev) {
                return ev.weight == 0;
            }),
        ev_collisions.end());
    // Same for edge-edge collisions.
    ee_collisions.erase(
        std::remove_if(
            ee_collisions.begin(), ee_collisions.end(),
            [&](const EdgeEdgeNormalCollision& ee) { return ee.weight == 0; }),
        ee_collisions.end());
}

} // namespace ipc