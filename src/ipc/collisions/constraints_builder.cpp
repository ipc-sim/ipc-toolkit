#include "constraints_builder.hpp"

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/point_triangle.hpp>

namespace ipc {

CollisionConstraintsBuilder::CollisionConstraintsBuilder(
    const bool use_convergent_formulation, const bool compute_shape_derivatives)
    : use_convergent_formulation(use_convergent_formulation)
    , compute_shape_derivatives(compute_shape_derivatives)
{
}

// ============================================================================

void CollisionConstraintsBuilder::add_edge_vertex_constraints(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const std::vector<EdgeVertexCandidate>& candidates,
    const std::function<bool(double)>& is_active,
    const size_t start_i,
    const size_t end_i)
{
    auto& C_ev = constraints.ev_constraints;
    const Eigen::MatrixXi& E = mesh.edges();

    for (size_t i = start_i; i < end_i; i++) {
        const auto& [ei, vi] = candidates[i];
        long e0i = E(ei, 0), e1i = E(ei, 1);

        PointEdgeDistanceType dtype =
            point_edge_distance_type(V.row(vi), V.row(e0i), V.row(e1i));

        double distance_sqr =
            point_edge_distance(V.row(vi), V.row(e0i), V.row(e1i), dtype);

        if (!is_active(distance_sqr))
            continue;

        // ÷ 2 to handle double counting for correct integration
        const double weight =
            use_convergent_formulation ? (mesh.vertex_area(vi) / 2) : 1;

        Eigen::SparseVector<double> weight_gradient;
        if (constraints.compute_shape_derivatives) {
            weight_gradient = constraints.use_convergent_formulation
                ? (mesh.vertex_area_gradient(vi) / 2)
                : Eigen::SparseVector<double>(V.size());
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
            C_ev.emplace_back(ei, vi);
            C_ev.back().weight = weight;
            C_ev.back().weight_gradient = weight_gradient;
            ev_to_index.emplace(C_ev.back(), C_ev.size() - 1);
            break;
        }
    }
}

void CollisionConstraintsBuilder::add_edge_edge_constraints(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const std::vector<EdgeEdgeCandidate>& candidates,
    const std::function<bool(double)>& is_active,
    const size_t start_i,
    const size_t end_i)
{
    auto& C_ee = constraints.ee_constraints;
    const Eigen::MatrixXd& V_rest = mesh.vertices_at_rest();
    const Eigen::MatrixXi& E = mesh.edges();

    for (size_t i = start_i; i < end_i; i++) {
        const auto& [eai, ebi] = candidates[i];
        long ea0i = E(eai, 0), ea1i = E(eai, 1);
        long eb0i = E(ebi, 0), eb1i = E(ebi, 1);

        EdgeEdgeDistanceType dtype = edge_edge_distance_type(
            V.row(ea0i), V.row(ea1i), V.row(eb0i), V.row(eb1i));

        double distance_sqr = edge_edge_distance(
            V.row(ea0i), V.row(ea1i), V.row(eb0i), V.row(eb1i), dtype);

        if (!is_active(distance_sqr))
            continue;

        // ÷ 4 to handle double counting and PT + EE for correct integration.
        // Sum edge areas because duplicate edge candidates were removed.
        const double weight = use_convergent_formulation
            ? ((mesh.edge_area(eai) + mesh.edge_area(ebi)) / 4)
            : 1;

        Eigen::SparseVector<double> weight_gradient;
        if (compute_shape_derivatives) {
            weight_gradient = use_convergent_formulation
                ? ((mesh.edge_area_gradient(eai) + mesh.edge_area_gradient(ebi))
                   / 4)
                : Eigen::SparseVector<double>(V.size());
        }

        double eps_x = edge_edge_mollifier_threshold(
            V_rest.row(ea0i), V_rest.row(ea1i), //
            V_rest.row(eb0i), V_rest.row(eb1i));
        double ee_cross_norm_sqr = edge_edge_cross_squarednorm(
            V.row(ea0i), V.row(ea1i), V.row(eb0i), V.row(eb1i));
        if (ee_cross_norm_sqr < eps_x) {
            // NOTE: This may not actually be the distance type, but all EE
            // pairs requiring mollification must be mollified later.
            dtype = EdgeEdgeDistanceType::EA_EB;
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
            C_ee.emplace_back(eai, ebi, eps_x);
            C_ee.back().weight = weight;
            C_ee.back().weight_gradient = weight_gradient;
            break;
        }
    }
}

void CollisionConstraintsBuilder::add_face_vertex_constraints(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const std::vector<FaceVertexCandidate>& candidates,
    const std::function<bool(double)>& is_active,
    const size_t start_i,
    const size_t end_i)
{
    auto& C_fv = constraints.fv_constraints;
    const Eigen::MatrixXi& F = mesh.faces();
    const Eigen::MatrixXi& F2E = mesh.faces_to_edges();

    for (size_t i = start_i; i < end_i; i++) {
        const auto& [fi, vi] = candidates[i];
        long f0i = F(fi, 0), f1i = F(fi, 1), f2i = F(fi, 2);

        // Compute distance type
        PointTriangleDistanceType dtype = point_triangle_distance_type(
            V.row(vi), V.row(f0i), V.row(f1i), V.row(f2i));

        double distance_sqr = point_triangle_distance(
            V.row(vi), V.row(f0i), V.row(f1i), V.row(f2i), dtype);

        if (!is_active(distance_sqr))
            continue;

        // ÷ 4 to handle double counting and PT + EE for correct integration
        const double weight =
            use_convergent_formulation ? (mesh.vertex_area(vi) / 4) : 1;

        Eigen::SparseVector<double> weight_gradient;
        if (compute_shape_derivatives) {
            weight_gradient = use_convergent_formulation
                ? (mesh.vertex_area_gradient(vi) / 4)
                : Eigen::SparseVector<double>(V.size());
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
            add_edge_vertex_constraint(F2E(fi, 0), vi, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::P_E1:
            add_edge_vertex_constraint(F2E(fi, 1), vi, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::P_E2:
            add_edge_vertex_constraint(F2E(fi, 2), vi, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::P_T:
            C_fv.emplace_back(fi, vi);
            C_fv.back().weight = weight;
            C_fv.back().weight_gradient = weight_gradient;
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
    unordered_map<VertexVertexConstraint, long>& vv_to_index,
    std::vector<VertexVertexConstraint>& vv_constraints)
{
    VertexVertexConstraint vv_constraint(v0i, v1i);
    auto found_item = vv_to_index.find(vv_constraint);
    if (found_item != vv_to_index.end()) {
        // Constraint already exists, so increase weight
        vv_constraints[found_item->second].weight += weight;
        vv_constraints[found_item->second].weight_gradient += weight_gradient;
    } else {
        // New constraint, so add it to the end of vv_constraints
        vv_to_index.emplace(vv_constraint, vv_constraints.size());
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
    unordered_map<EdgeVertexConstraint, long>& ev_to_index,
    std::vector<EdgeVertexConstraint>& ev_constraints)
{
    EdgeVertexConstraint ev_constraint(ei, vi);
    auto found_item = ev_to_index.find(ev_constraint);
    if (found_item != ev_to_index.end()) {
        // Constraint already exists, so increase weight
        ev_constraints[found_item->second].weight += weight;
        ev_constraints[found_item->second].weight_gradient += weight_gradient;
    } else {
        // New constraint, so add it to the end of vv_constraints
        ev_to_index.emplace(ev_constraint, ev_constraints.size());
        ev_constraints.push_back(ev_constraint);
        ev_constraints.back().weight = weight;
        ev_constraints.back().weight_gradient = weight_gradient;
    }
}

// ============================================================================

void CollisionConstraintsBuilder::merge(
    const tbb::enumerable_thread_specific<CollisionConstraintsBuilder>&
        local_storage,
    CollisionConstraints& constraints)
{
    unordered_map<VertexVertexConstraint, long> vv_to_index;
    unordered_map<EdgeVertexConstraint, long> ev_to_index;
    auto& vv_constraints = constraints.vv_constraints;
    auto& ev_constraints = constraints.ev_constraints;
    auto& ee_constraints = constraints.ee_constraints;
    auto& fv_constraints = constraints.fv_constraints;

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
            vv_to_index = builder.vv_to_index;
            vv_constraints.insert(
                vv_constraints.end(), local_constraints.vv_constraints.begin(),
                local_constraints.vv_constraints.end());
        } else {
            for (const auto& vv : local_constraints.vv_constraints) {
                add_vertex_vertex_constraint(
                    vv.vertex0_index, vv.vertex1_index, vv.weight,
                    vv.weight_gradient, vv_to_index, vv_constraints);
            }
        }

        if (ev_constraints.empty()) {
            ev_to_index = builder.ev_to_index;
            ev_constraints.insert(
                ev_constraints.end(), local_constraints.ev_constraints.begin(),
                local_constraints.ev_constraints.end());
        } else {
            for (const auto& ev : local_constraints.ev_constraints) {
                add_edge_vertex_constraint(
                    ev.edge_index, ev.vertex_index, ev.weight,
                    ev.weight_gradient, ev_to_index, ev_constraints);
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