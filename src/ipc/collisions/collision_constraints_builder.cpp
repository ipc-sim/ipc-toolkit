#include "collision_constraints_builder.hpp"

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/point_triangle.hpp>

namespace ipc {

CollisionConstraintsBuilder::CollisionConstraintsBuilder(
    const bool _use_convergent_formulation,
    const bool _should_compute_weight_gradient)
    : use_convergent_formulation(_use_convergent_formulation)
    , should_compute_weight_gradient(_should_compute_weight_gradient)
{
}

// ============================================================================

void CollisionConstraintsBuilder::add_vertex_vertex_constraints(
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
        const double weight = use_convergent_formulation
            ? ((mesh.vertex_area(vi) + mesh.vertex_area(vj)) / 2)
            : 1;

        Eigen::SparseVector<double> weight_gradient;
        if (should_compute_weight_gradient) {
            weight_gradient = use_convergent_formulation
                ? ((mesh.vertex_area_gradient(vi)
                    + mesh.vertex_area_gradient(vj))
                   / 2)
                : Eigen::SparseVector<double>(vertices.size());
        }

        VertexVertexConstraint vv_constraint(vi, vj, weight, weight_gradient);
        vv_to_id.emplace(vv_constraint, vv_constraints.size());
        vv_constraints.push_back(vv_constraint);
    }
}

void CollisionConstraintsBuilder::add_edge_vertex_constraints(
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
            use_convergent_formulation ? (mesh.vertex_area(vi) / 2) : 1;

        Eigen::SparseVector<double> weight_gradient;
        if (should_compute_weight_gradient) {
            weight_gradient = use_convergent_formulation
                ? (mesh.vertex_area_gradient(vi) / 2)
                : Eigen::SparseVector<double>(vertices.size());
        }

        add_edge_vertex_constraint(
            mesh, candidates[i], dtype, weight, weight_gradient);
    }
}

void CollisionConstraintsBuilder::add_edge_vertex_constraint(
    const CollisionMesh& mesh,
    const EdgeVertexCandidate& candidate,
    const PointEdgeDistanceType dtype,
    const double weight,
    const Eigen::SparseVector<double>& weight_gradient)
{
    const auto& [ei, vi] = candidate;

    switch (dtype) {
    case PointEdgeDistanceType::P_E0:
        add_vertex_vertex_constraint(
            vi, mesh.edges()(ei, 0), weight, weight_gradient);
        break;

    case PointEdgeDistanceType::P_E1:
        add_vertex_vertex_constraint(
            vi, mesh.edges()(ei, 1), weight, weight_gradient);
        break;

    case PointEdgeDistanceType::P_E:
        add_edge_vertex_constraint(ei, vi, weight, weight_gradient);
        break;

    case PointEdgeDistanceType::AUTO:
    default:
        assert(false);
        break;
    }
}

void CollisionConstraintsBuilder::add_edge_edge_constraints(
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
        const double weight = use_convergent_formulation
            ? ((mesh.edge_area(eai) + mesh.edge_area(ebi)) / 4)
            : 1;

        Eigen::SparseVector<double> weight_gradient;
        if (should_compute_weight_gradient) {
            weight_gradient = use_convergent_formulation
                ? ((mesh.edge_area_gradient(eai) + mesh.edge_area_gradient(ebi))
                   / 4)
                : Eigen::SparseVector<double>(vertices.size());
        }

        switch (dtype) {
        case EdgeEdgeDistanceType::EA0_EB0:
            add_vertex_vertex_constraint(ea0i, eb0i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA0_EB1:
            add_vertex_vertex_constraint(ea0i, eb1i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA1_EB0:
            add_vertex_vertex_constraint(ea1i, eb0i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA1_EB1:
            add_vertex_vertex_constraint(ea1i, eb1i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA_EB0:
            add_edge_vertex_constraint(eai, eb0i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA_EB1:
            add_edge_vertex_constraint(eai, eb1i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA0_EB:
            add_edge_vertex_constraint(ebi, ea0i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA1_EB:
            add_edge_vertex_constraint(ebi, ea1i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA_EB:
            ee_constraints.emplace_back(
                eai, ebi, eps_x, weight, weight_gradient, actual_dtype);
            ee_to_id.emplace(ee_constraints.back(), ee_constraints.size() - 1);
            break;

        case EdgeEdgeDistanceType::AUTO:
        default:
            assert(false);
            break;
        }
    }
}

void CollisionConstraintsBuilder::add_face_vertex_constraints(
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
            use_convergent_formulation ? (mesh.vertex_area(vi) / 4) : 1;

        Eigen::SparseVector<double> weight_gradient;
        if (should_compute_weight_gradient) {
            weight_gradient = use_convergent_formulation
                ? (mesh.vertex_area_gradient(vi) / 4)
                : Eigen::SparseVector<double>(vertices.size());
        }

        switch (dtype) {
        case PointTriangleDistanceType::P_T0:
            add_vertex_vertex_constraint(vi, f0i, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::P_T1:
            add_vertex_vertex_constraint(vi, f1i, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::P_T2:
            add_vertex_vertex_constraint(vi, f2i, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::P_E0:
            add_edge_vertex_constraint(
                mesh.faces_to_edges()(fi, 0), vi, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::P_E1:
            add_edge_vertex_constraint(
                mesh.faces_to_edges()(fi, 1), vi, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::P_E2:
            add_edge_vertex_constraint(
                mesh.faces_to_edges()(fi, 2), vi, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::P_T:
            fv_constraints.emplace_back(fi, vi, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::AUTO:
        default:
            assert(false);
            break;
        }
    }
}

// ============================================================================

void CollisionConstraintsBuilder::
    add_edge_vertex_negative_vertex_vertex_constraints(
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
                * (use_convergent_formulation ? (mesh.vertex_area(vi) / 2) : 1);

            if (should_compute_weight_gradient && use_convergent_formulation) {
                weight_gradient += (1 - incident_edge_amt) / 2.0
                    * mesh.vertex_area_gradient(vi);
            }
        }
    };

    for (size_t i = start_i; i < end_i; i++) {
        const auto& [vi, vj] = candidates[i];
        assert(vi != vj);

        double weight = 0;
        Eigen::SparseVector<double> weight_gradient;
        if (should_compute_weight_gradient) {
            weight_gradient = Eigen::SparseVector<double>(vertices.size());
        }

        add_weight(vi, vj, weight, weight_gradient);
        add_weight(vj, vi, weight, weight_gradient);

        if (weight != 0) {
            add_vertex_vertex_constraint(vi, vj, weight, weight_gradient);
        }
    }
}

void CollisionConstraintsBuilder::
    add_face_vertex_positive_vertex_vertex_constraints(
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
        weight += use_convergent_formulation ? (mesh.vertex_area(vi) / 4) : 1;

        if (should_compute_weight_gradient && use_convergent_formulation) {
            weight_gradient += mesh.vertex_area_gradient(vi) / 4;
        }
    };

    for (size_t i = start_i; i < end_i; i++) {
        const auto& [vi, vj] = candidates[i];
        assert(vi != vj);

        double weight = 0;
        Eigen::SparseVector<double> weight_gradient;
        if (should_compute_weight_gradient) {
            weight_gradient = Eigen::SparseVector<double>(vertices.size());
        }

        add_weight(vi, vj, weight, weight_gradient);
        add_weight(vj, vi, weight, weight_gradient);

        if (weight != 0) {
            add_vertex_vertex_constraint(vi, vj, weight, weight_gradient);
        }
    }
}

void CollisionConstraintsBuilder::
    add_face_vertex_negative_edge_vertex_constraints(
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
                * (use_convergent_formulation ? (mesh.vertex_area(vi) / 4) : 1);

            Eigen::SparseVector<double> weight_gradient;
            if (should_compute_weight_gradient) {
                weight_gradient = use_convergent_formulation
                    ? ((1 - incident_triangle_amt) / 4.0
                       * mesh.vertex_area_gradient(vi))
                    : Eigen::SparseVector<double>(vertices.size());
            }

            add_edge_vertex_constraint(
                mesh, candidates[i],
                point_edge_distance_type(
                    vertices.row(vi), vertices.row(mesh.edges()(ei, 0)),
                    vertices.row(mesh.edges()(ei, 1))),
                weight, weight_gradient);
        }
    }
}

void CollisionConstraintsBuilder::
    add_edge_edge_negative_edge_vertex_constraints(
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
            use_convergent_formulation ? (-0.25 * mesh.edge_area(ea)) : -1;
        Eigen::SparseVector<double> weight_gradient;
        if (should_compute_weight_gradient) {
            weight_gradient = use_convergent_formulation
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

            // Add mollified EE constraint with specified distance type
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

            add_edge_edge_constraint(
                ea, eb, eps_x, weight, weight_gradient, ee_dtype);
        }

        if (nonmollified_incident_edge_amt == 1) {
            continue; // no constraint to add because (ρ(x) - 1) = 0
        }
        // if nonmollified_incident_edge_amt == 0, then we need to explicitly
        // add a positive constraint.
        add_edge_vertex_constraint(
            mesh, candidates[i], dtype,
            (nonmollified_incident_edge_amt - 1) * weight,
            (nonmollified_incident_edge_amt - 1) * weight_gradient);
    }
}

// ============================================================================

void CollisionConstraintsBuilder::add_vertex_vertex_constraint(
    const VertexVertexConstraint& vv_constraint,
    unordered_map<VertexVertexConstraint, long>& vv_to_id,
    std::vector<VertexVertexConstraint>& vv_constraints)
{
    auto found_item = vv_to_id.find(vv_constraint);
    if (found_item != vv_to_id.end()) {
        // Constraint already exists, so increase weight
        vv_constraints[found_item->second].weight += vv_constraint.weight;
        vv_constraints[found_item->second].weight_gradient +=
            vv_constraint.weight_gradient;
    } else {
        // New constraint, so add it to the end of vv_constraints
        vv_to_id.emplace(vv_constraint, vv_constraints.size());
        vv_constraints.push_back(vv_constraint);
    }
}

void CollisionConstraintsBuilder::add_edge_vertex_constraint(
    const EdgeVertexConstraint& ev_constraint,
    unordered_map<EdgeVertexConstraint, long>& ev_to_id,
    std::vector<EdgeVertexConstraint>& ev_constraints)
{
    auto found_item = ev_to_id.find(ev_constraint);
    if (found_item != ev_to_id.end()) {
        // Constraint already exists, so increase weight
        ev_constraints[found_item->second].weight += ev_constraint.weight;
        ev_constraints[found_item->second].weight_gradient +=
            ev_constraint.weight_gradient;
    } else {
        // New constraint, so add it to the end of ev_constraints
        ev_to_id.emplace(ev_constraint, ev_constraints.size());
        ev_constraints.push_back(ev_constraint);
    }
}

void CollisionConstraintsBuilder::add_edge_edge_constraint(
    const EdgeEdgeConstraint& ee_constraint,
    unordered_map<EdgeEdgeConstraint, long>& ee_to_id,
    std::vector<EdgeEdgeConstraint>& ee_constraints)
{
    auto found_item = ee_to_id.find(ee_constraint);
    if (found_item != ee_to_id.end()) {
        // Constraint already exists, so increase weight
        assert(ee_constraint == ee_constraints[found_item->second]);
        ee_constraints[found_item->second].weight += ee_constraint.weight;
        ee_constraints[found_item->second].weight_gradient +=
            ee_constraint.weight_gradient;
    } else {
        // New constraint, so add it to the end of ee_constraints
        ee_to_id.emplace(ee_constraint, ee_constraints.size());
        ee_constraints.push_back(ee_constraint);
    }
}

// ============================================================================

void CollisionConstraintsBuilder::merge(
    const tbb::enumerable_thread_specific<CollisionConstraintsBuilder>&
        local_storage,
    CollisionConstraints& merged_constraints)
{
    unordered_map<VertexVertexConstraint, long> vv_to_id;
    unordered_map<EdgeVertexConstraint, long> ev_to_id;
    unordered_map<EdgeEdgeConstraint, long> ee_to_id;
    auto& vv_constraints = merged_constraints.vv_constraints;
    auto& ev_constraints = merged_constraints.ev_constraints;
    auto& ee_constraints = merged_constraints.ee_constraints;
    auto& fv_constraints = merged_constraints.fv_constraints;

    // size up the hash items
    size_t n_vv = 0, n_ev = 0, n_ee = 0, n_fv = 0;
    for (const auto& storage : local_storage) {
        // This is an conservative estimate
        n_vv += storage.vv_constraints.size();
        n_ev += storage.ev_constraints.size();
        n_ee += storage.ee_constraints.size();
        n_fv += storage.fv_constraints.size();
    }
    vv_constraints.reserve(n_vv);
    ev_constraints.reserve(n_ev);
    ee_constraints.reserve(n_ee);
    fv_constraints.reserve(n_fv);

    // merge
    for (const auto& builder : local_storage) {
        if (vv_constraints.empty()) {
            vv_to_id = builder.vv_to_id;
            vv_constraints = builder.vv_constraints;
        } else {
            for (const auto& vv : builder.vv_constraints) {
                add_vertex_vertex_constraint(vv, vv_to_id, vv_constraints);
            }
        }

        if (ev_constraints.empty()) {
            ev_to_id = builder.ev_to_id;
            ev_constraints = builder.ev_constraints;
        } else {
            for (const auto& ev : builder.ev_constraints) {
                add_edge_vertex_constraint(ev, ev_to_id, ev_constraints);
            }
        }

        if (ee_constraints.empty()) {
            ee_to_id = builder.ee_to_id;
            ee_constraints = builder.ee_constraints;
        } else {
            for (const auto& ee : builder.ee_constraints) {
                add_edge_edge_constraint(ee, ee_to_id, ee_constraints);
            }
        }

        fv_constraints.insert(
            fv_constraints.end(), builder.fv_constraints.begin(),
            builder.fv_constraints.end());
    }

    // If positive and negative vertex-vertex constraints cancel out, remove
    // them. This can happen when edge-vertex constraints reduce to
    // vertex-vertex constraints. This will avoid unnecessary computation.
    vv_constraints.erase(
        std::remove_if(
            vv_constraints.begin(), vv_constraints.end(),
            [&](const VertexVertexConstraint& vv) { return vv.weight == 0; }),
        vv_constraints.end());
    // Same for edge-vertex constraints.
    ev_constraints.erase(
        std::remove_if(
            ev_constraints.begin(), ev_constraints.end(),
            [&](const EdgeVertexConstraint& ev) { return ev.weight == 0; }),
        ev_constraints.end());
    // Same for edge-edge constraints.
    ee_constraints.erase(
        std::remove_if(
            ee_constraints.begin(), ee_constraints.end(),
            [&](const EdgeEdgeConstraint& ee) { return ee.weight == 0; }),
        ee_constraints.end());
}

} // namespace ipc