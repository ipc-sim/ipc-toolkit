#include "collision_constraints_builder.hpp"

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/point_triangle.hpp>

namespace ipc {

CollisionConstraintsBuilder::CollisionConstraintsBuilder(
    const CollisionConstraints& empty_constraints)
{
    assert(empty_constraints.empty());
    constraints = empty_constraints;
}

// ============================================================================

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
        const long e0i = mesh.edges()(ei, 0), e1i = mesh.edges()(ei, 1);

        const auto [v, e0, e1, _] =
            candidates[i].vertices(vertices, mesh.edges(), mesh.faces());
        const PointEdgeDistanceType dtype = point_edge_distance_type(v, e0, e1);
        const double distance_sqr = point_edge_distance(v, e0, e1, dtype);

        if (!is_active(distance_sqr))
            continue;

        // ÷ 2 to handle double counting for correct integration
        const double weight =
            use_convergent_formulation() ? (mesh.vertex_area(vi) / 2) : 1;

        Eigen::SparseVector<double> weight_gradient;
        if (should_compute_weight_gradient()) {
            weight_gradient = use_convergent_formulation()
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

        EdgeEdgeDistanceType dtype =
            edge_edge_distance_type(ea0, ea1, eb0, eb1);

        const double distance_sqr =
            edge_edge_distance(ea0, ea1, eb0, eb1, dtype);

        if (!is_active(distance_sqr))
            continue;

        const double eps_x = edge_edge_mollifier_threshold(
            mesh.rest_positions().row(ea0i), mesh.rest_positions().row(ea1i),
            mesh.rest_positions().row(eb0i), mesh.rest_positions().row(eb1i));

        const double ee_cross_norm_sqr =
            edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1);

        if (ee_cross_norm_sqr < eps_x) {
            // NOTE: This may not actually be the distance type, but all EE
            // pairs requiring mollification must be mollified later.
            dtype = EdgeEdgeDistanceType::EA_EB;
        }

        // ÷ 4 to handle double counting and PT + EE for correct integration.
        // Sum edge areas because duplicate edge candidates were removed.
        const double weight = use_convergent_formulation()
            ? ((mesh.edge_area(eai) + mesh.edge_area(ebi)) / 4)
            : 1;

        Eigen::SparseVector<double> weight_gradient;
        if (should_compute_weight_gradient()) {
            weight_gradient = use_convergent_formulation()
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
            constraints.ee_constraints.emplace_back(eai, ebi, eps_x);
            constraints.ee_constraints.back().weight = weight;
            constraints.ee_constraints.back().weight_gradient = weight_gradient;
            break;

        case EdgeEdgeDistanceType::AUTO:
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
            use_convergent_formulation() ? (mesh.vertex_area(vi) / 4) : 1;

        Eigen::SparseVector<double> weight_gradient;
        if (should_compute_weight_gradient()) {
            weight_gradient = use_convergent_formulation()
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
            constraints.fv_constraints.emplace_back(fi, vi);
            constraints.fv_constraints.back().weight = weight;
            constraints.fv_constraints.back().weight_gradient = weight_gradient;
            break;

        case PointTriangleDistanceType::AUTO:
            assert(false);
            break;
        }
    }
}

// ============================================================================

void CollisionConstraintsBuilder::add_vertex_vertex_negative_constraints(
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
                * (use_convergent_formulation() ? (mesh.vertex_area(vi) / 2)
                                                : 1);

            if (should_compute_weight_gradient()
                && use_convergent_formulation()) {
                weight_gradient += (1 - incident_edge_amt)
                    * (mesh.vertex_area_gradient(vi) / 2);
            }
        }
    };

    for (size_t i = start_i; i < end_i; i++) {
        const auto& [vi, vj] = candidates[i];
        assert(vi != vj);

        double weight = 0;
        Eigen::SparseVector<double> weight_gradient;
        if (should_compute_weight_gradient()) {
            weight_gradient = Eigen::SparseVector<double>(vertices.size());
        }

        add_weight(vi, vj, weight, weight_gradient);
        add_weight(vj, vi, weight, weight_gradient);

        if (weight != 0) {
            add_vertex_vertex_constraint(vi, vj, weight, weight_gradient);
        }
    }
}

void CollisionConstraintsBuilder::add_vertex_vertex_positive_constraints(
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
        weight += use_convergent_formulation() ? (mesh.vertex_area(vi) / 4) : 1;

        if (should_compute_weight_gradient() && use_convergent_formulation()) {
            weight_gradient += mesh.vertex_area_gradient(vi) / 4;
        }
    };

    for (size_t i = start_i; i < end_i; i++) {
        const auto& [vi, vj] = candidates[i];
        assert(vi != vj);

        double weight = 0;
        Eigen::SparseVector<double> weight_gradient;
        if (should_compute_weight_gradient()) {
            weight_gradient = Eigen::SparseVector<double>(vertices.size());
        }

        add_weight(vi, vj, weight, weight_gradient);
        add_weight(vj, vi, weight, weight_gradient);

        if (weight != 0) {
            add_vertex_vertex_constraint(vi, vj, weight, weight_gradient);
        }
    }
}

void CollisionConstraintsBuilder::add_edge_vertex_negative_constraints(
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
                * (use_convergent_formulation() ? (mesh.vertex_area(vi) / 4)
                                                : 1);

            Eigen::SparseVector<double> weight_gradient;
            if (should_compute_weight_gradient()
                && use_convergent_formulation()) {
                weight_gradient = (1 - incident_triangle_amt)
                    * (mesh.vertex_area_gradient(vi) / 4);
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

void CollisionConstraintsBuilder::add_edge_vertex_negative_constraints(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<
        std::tuple<EdgeVertexCandidate, double, Eigen::SparseVector<double>>>&
        candidates,
    const size_t start_i,
    const size_t end_i)
{
    for (size_t i = start_i; i < end_i; i++) {
        const auto& [ei, vi] = std::get<0>(candidates[i]);
        const int e0i = mesh.edges()(ei, 0), e1i = mesh.edges()(ei, 1);
        assert(vi != e0i && vi != e1i);

        // TODO: distinguish mollified vs non-mollified
        const auto& incident_vertices = mesh.vertex_vertex_adjacencies()[vi];
        const int incident_edge_amt = incident_vertices.size()
            - int(incident_vertices.find(mesh.edges()(e0i, 0))
                  != incident_vertices.end())
            - int(incident_vertices.find(mesh.edges()(e1i, 1))
                  != incident_vertices.end());

        if (incident_edge_amt > 1) {
            // ÷ 4 to handle double counting and PT + EE for correct integration
            const double weight = (1 - incident_edge_amt)
                * (use_convergent_formulation()
                       ? ((mesh.edge_area(ei) + std::get<1>(candidates[i])) / 4)
                       : 1);

            Eigen::SparseVector<double> weight_gradient;
            if (should_compute_weight_gradient()
                && use_convergent_formulation()) {
                weight_gradient = (1 - incident_edge_amt)
                    * (mesh.edge_area_gradient(ei) + std::get<2>(candidates[i]))
                    / 4;
            }

            add_edge_vertex_constraint(
                mesh, std::get<0>(candidates[i]),
                point_edge_distance_type(
                    vertices.row(vi), vertices.row(mesh.edges()(ei, 0)),
                    vertices.row(mesh.edges()(ei, 1))),
                weight, weight_gradient);
        }
    }
}

// ============================================================================

void CollisionConstraintsBuilder::add_vertex_vertex_constraint(
    const long v0i,
    const long v1i,
    const double weight,
    const Eigen::SparseVector<double>& weight_gradient,
    unordered_map<VertexVertexConstraint, long>& vv_to_id,
    std::vector<VertexVertexConstraint>& vv_constraints)
{
    VertexVertexConstraint vv_constraint(v0i, v1i);
    auto found_item = vv_to_id.find(vv_constraint);
    if (found_item != vv_to_id.end()) {
        // Constraint already exists, so increase weight
        vv_constraints[found_item->second].weight += weight;
        vv_constraints[found_item->second].weight_gradient += weight_gradient;
    } else {
        // New constraint, so add it to the end of vv_constraints
        vv_to_id.emplace(vv_constraint, vv_constraints.size());
        vv_constraints.push_back(vv_constraint);
        vv_constraints.back().weight = weight;
        vv_constraints.back().weight_gradient = weight_gradient;
    }
}

void CollisionConstraintsBuilder::add_edge_vertex_constraint(
    const long ei,
    const long vi,
    const double weight,
    const Eigen::SparseVector<double>& weight_gradient,
    unordered_map<EdgeVertexConstraint, long>& ev_to_id,
    std::vector<EdgeVertexConstraint>& ev_constraints)
{
    EdgeVertexConstraint ev_constraint(ei, vi);
    auto found_item = ev_to_id.find(ev_constraint);
    if (found_item != ev_to_id.end()) {
        // Constraint already exists, so increase weight
        ev_constraints[found_item->second].weight += weight;
        ev_constraints[found_item->second].weight_gradient += weight_gradient;
    } else {
        // New constraint, so add it to the end of ev_constraints
        ev_to_id.emplace(ev_constraint, ev_constraints.size());
        ev_constraints.push_back(ev_constraint);
        ev_constraints.back().weight = weight;
        ev_constraints.back().weight_gradient = weight_gradient;
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
    auto& vv_constraints = merged_constraints.vv_constraints;
    auto& ev_constraints = merged_constraints.ev_constraints;
    auto& ee_constraints = merged_constraints.ee_constraints;
    auto& fv_constraints = merged_constraints.fv_constraints;

    // size up the hash items
    size_t n_vv = 0, n_ev = 0, n_ee = 0, n_fv = 0;
    for (const auto& storage : local_storage) {
        // This is an conservative estimate
        n_vv += storage.constraints.vv_constraints.size();
        n_ev += storage.constraints.ev_constraints.size();
        n_ee += storage.constraints.ee_constraints.size();
        n_fv += storage.constraints.fv_constraints.size();
    }
    vv_constraints.reserve(n_vv);
    ev_constraints.reserve(n_ev);
    ee_constraints.reserve(n_ee);
    fv_constraints.reserve(n_fv);

    // merge
    for (const auto& builder : local_storage) {
        const auto& local_constraints = builder.constraints;

        if (vv_constraints.empty()) {
            vv_to_id = builder.vv_to_id;
            vv_constraints.insert(
                vv_constraints.end(), local_constraints.vv_constraints.begin(),
                local_constraints.vv_constraints.end());
        } else {
            for (const auto& vv : local_constraints.vv_constraints) {
                add_vertex_vertex_constraint(
                    vv.vertex0_id, vv.vertex1_id, vv.weight, vv.weight_gradient,
                    vv_to_id, vv_constraints);
            }
        }

        if (ev_constraints.empty()) {
            ev_to_id = builder.ev_to_id;
            ev_constraints.insert(
                ev_constraints.end(), local_constraints.ev_constraints.begin(),
                local_constraints.ev_constraints.end());
        } else {
            for (const auto& ev : local_constraints.ev_constraints) {
                add_edge_vertex_constraint(
                    ev.edge_id, ev.vertex_id, ev.weight, ev.weight_gradient,
                    ev_to_id, ev_constraints);
            }
        }

        ee_constraints.insert(
            ee_constraints.end(), local_constraints.ee_constraints.begin(),
            local_constraints.ee_constraints.end());
        fv_constraints.insert(
            fv_constraints.end(), local_constraints.fv_constraints.begin(),
            local_constraints.fv_constraints.end());
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
}

} // namespace ipc