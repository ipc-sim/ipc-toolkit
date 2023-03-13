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
        PointEdgeDistanceType dtype = point_edge_distance_type(v, e0, e1);
        double distance_sqr = point_edge_distance(v, e0, e1, dtype);

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

        switch (dtype) {
        case PointEdgeDistanceType::P_E0:
            add_vertex_vertex_constraint(vi, e0i, weight, weight_gradient);
            break;

        case PointEdgeDistanceType::P_E1:
            add_vertex_vertex_constraint(vi, e1i, weight, weight_gradient);
            break;

        case PointEdgeDistanceType::P_E:
            // ev_candidates is a set, so no duplicate EV CollisionConstraints
            constraints.ev_constraints.emplace_back(ei, vi);
            constraints.ev_constraints.back().weight = weight;
            constraints.ev_constraints.back().weight_gradient = weight_gradient;
            ev_to_id.emplace(
                constraints.ev_constraints.back(),
                constraints.ev_constraints.size() - 1);
            break;

        case PointEdgeDistanceType::AUTO:
            assert(false);
            break;
        }
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
        // New constraint, so add it to the end of vv_constraints
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
}

} // namespace ipc