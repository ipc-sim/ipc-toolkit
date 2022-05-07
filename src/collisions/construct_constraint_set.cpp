#include <ipc/collisions/construct_constraint_set.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/utils/unordered_map_and_set.hpp>

namespace ipc {

void construct_constraint_set(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const double dhat,
    Constraints& constraint_set,
    const double dmin,
    const BroadPhaseMethod method)
{
    assert(V.rows() == mesh.num_vertices());

    double inflation_radius = (dhat + dmin) / 1.99; // Conservative inflation

    Candidates candidates;
    construct_collision_candidates(
        mesh, V, candidates, inflation_radius, method);

    construct_constraint_set(candidates, mesh, V, dhat, constraint_set, dmin);
}

void construct_constraint_set(
    const Candidates& candidates,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const double dhat,
    Constraints& constraint_set,
    const double dmin)
{
    assert(V.rows() == mesh.num_vertices());

    constraint_set.clear();

    // Cull the candidates by measuring the distance and dropping those that are
    // greater than dhat.
    const double offset_sqr = (dmin + dhat) * (dmin + dhat);
    auto is_active = [&](double distance_sqr) {
        return distance_sqr < offset_sqr;
    };

    tbb::enumerable_thread_specific<ConstraintBuilder> storage;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), candidates.ev_candidates.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            edge_vertex_candiates_to_constraints(
                mesh, V, candidates.ev_candidates, is_active, r.begin(),
                r.end(), storage.local());
        });

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), candidates.ee_candidates.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            edge_edge_candiates_to_constraints(
                mesh, V, candidates.ee_candidates, is_active, r.begin(),
                r.end(), storage.local());
        });

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), candidates.fv_candidates.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            face_vertex_candiates_to_constraints(
                mesh, V, candidates.fv_candidates, is_active, r.begin(),
                r.end(), storage.local());
        });

    merge_thread_local_constraints(storage, constraint_set);

    for (size_t ci = 0; ci < constraint_set.size(); ci++) {
        constraint_set[ci].minimum_distance = dmin;
    }
}

///////////////////////////////////////////////////////////////////////////////

template <typename Hash>
void add_vertex_vertex_constraint(
    std::vector<VertexVertexConstraint>& vv_constraints,
    unordered_map<VertexVertexConstraint, long, Hash>& vv_to_index,
    const long v0i,
    const long v1i,
    const double weight,
    const Eigen::SparseVector<double>& weight_gradient)
{
    VertexVertexConstraint vv_constraint(v0i, v1i);
    auto found_item = vv_to_index.find(vv_constraint);
    if (found_item != vv_to_index.end()) {
        // Constraint already exists, so increase weight
        vv_constraints[found_item->second].weight += weight;
#ifdef IPC_TOOLKIT_COMPUTE_SHAPE_DERIVATIVE
        vv_constraints[found_item->second].weight_gradient += weight_gradient;
#endif
    } else {
        // New constraint, so add it to the end of vv_constraints
        vv_to_index.emplace(vv_constraint, vv_constraints.size());
        vv_constraints.push_back(vv_constraint);
        vv_constraints.back().weight = weight;
#ifdef IPC_TOOLKIT_COMPUTE_SHAPE_DERIVATIVE
        vv_constraints.back().weight_gradient = weight_gradient;
#endif
    }
}

template <typename Hash>
void add_edge_vertex_constraint(
    std::vector<EdgeVertexConstraint>& ev_constraints,
    unordered_map<EdgeVertexConstraint, long, Hash>& ev_to_index,
    const long ei,
    const long vi,
    const double weight,
    const Eigen::SparseVector<double>& weight_gradient)
{
    EdgeVertexConstraint ev_constraint(ei, vi);
    auto found_item = ev_to_index.find(ev_constraint);
    if (found_item != ev_to_index.end()) {
        // Constraint already exists, so increase weight
        ev_constraints[found_item->second].weight += weight;
#ifdef IPC_TOOLKIT_COMPUTE_SHAPE_DERIVATIVE
        ev_constraints[found_item->second].weight_gradient += weight_gradient;
#endif
    } else {
        // New constraint, so add it to the end of vv_constraints
        ev_to_index.emplace(ev_constraint, ev_constraints.size());
        ev_constraints.push_back(ev_constraint);
        ev_constraints.back().weight = weight;
#ifdef IPC_TOOLKIT_COMPUTE_SHAPE_DERIVATIVE
        ev_constraints.back().weight_gradient = weight_gradient;
#endif
    }
}

void edge_vertex_candiates_to_constraints(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const std::vector<EdgeVertexCandidate>& candidates,
    const std::function<bool(double)>& is_active,
    const size_t start_i,
    const size_t end_i,
    ConstraintBuilder& constraint_builder)
{
    auto& [vv_to_index, ev_to_index, constraint_set] = constraint_builder;
    auto& [C_vv, C_ev, _, __, ___] = constraint_set;
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

#ifdef IPC_TOOLKIT_CONVERGENT
        // ÷ 2 to handle double counting for correct integration
        const double weight = mesh.point_area(vi) / 2;
        const Eigen::SparseVector<double> weight_gradient =
            mesh.point_area_gradient(vi) / 2;
#else
        const double weight = 1;
        const Eigen::SparseVector<double> weight_gradient(V.rows());
#endif

        switch (dtype) {
        case PointEdgeDistanceType::P_E0:
            add_vertex_vertex_constraint(
                C_vv, vv_to_index, vi, e0i, weight, weight_gradient);
            break;

        case PointEdgeDistanceType::P_E1:
            add_vertex_vertex_constraint(
                C_vv, vv_to_index, vi, e1i, weight, weight_gradient);
            break;

        case PointEdgeDistanceType::P_E:
            // ev_candidates is a set, so no duplicate EV
            // constraints
            C_ev.emplace_back(ei, vi);
            C_ev.back().weight = weight;
#ifdef IPC_TOOLKIT_COMPUTE_SHAPE_DERIVATIVE
            C_ev.back().weight_gradient = weight_gradient;
#endif
            ev_to_index.emplace(C_ev.back(), C_ev.size() - 1);
            break;
        }
    }
}

void edge_edge_candiates_to_constraints(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const std::vector<EdgeEdgeCandidate>& candidates,
    const std::function<bool(double)>& is_active,
    const size_t start_i,
    const size_t end_i,
    ConstraintBuilder& constraint_builder)
{
    auto& [vv_to_index, ev_to_index, constraint_set] = constraint_builder;
    auto& [C_vv, C_ev, C_ee, _, __] = constraint_set;
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

#ifdef IPC_TOOLKIT_CONVERGENT
        // ÷ 4 to handle double counting and PT + EE for correct integration.
        // Sum edge areas because duplicate edge candidates were removed.
        const double weight = (mesh.edge_area(eai) + mesh.edge_area(ebi)) / 4;
        const Eigen::SparseVector<double> weight_gradient =
            (mesh.edge_area_gradient(eai) + mesh.edge_area_gradient(ebi)) / 4;
#else
        const double weight = 1;
        const Eigen::SparseVector<double> weight_gradient(V.rows());
#endif

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
            add_vertex_vertex_constraint(
                C_vv, vv_to_index, ea0i, eb0i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA0_EB1:
            add_vertex_vertex_constraint(
                C_vv, vv_to_index, ea0i, eb1i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA1_EB0:
            add_vertex_vertex_constraint(
                C_vv, vv_to_index, ea1i, eb0i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA1_EB1:
            add_vertex_vertex_constraint(
                C_vv, vv_to_index, ea1i, eb1i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA_EB0:
            add_edge_vertex_constraint(
                C_ev, ev_to_index, eai, eb0i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA_EB1:
            add_edge_vertex_constraint(
                C_ev, ev_to_index, eai, eb1i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA0_EB:
            add_edge_vertex_constraint(
                C_ev, ev_to_index, ebi, ea0i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA1_EB:
            add_edge_vertex_constraint(
                C_ev, ev_to_index, ebi, ea1i, weight, weight_gradient);
            break;

        case EdgeEdgeDistanceType::EA_EB:
            C_ee.emplace_back(eai, ebi, eps_x);
            C_ee.back().weight = weight;
#ifdef IPC_TOOLKIT_COMPUTE_SHAPE_DERIVATIVE
            C_ee.back().weight_gradient = weight_gradient;
#endif
            break;
        }
    }
}

void face_vertex_candiates_to_constraints(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const std::vector<FaceVertexCandidate>& candidates,
    const std::function<bool(double)>& is_active,
    const size_t start_i,
    const size_t end_i,
    ConstraintBuilder& constraint_builder)
{
    auto& [vv_to_index, ev_to_index, constraint_set] = constraint_builder;
    auto& [C_vv, C_ev, _, C_fv, __] = constraint_set;
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

#ifdef IPC_TOOLKIT_CONVERGENT
        // ÷ 4 to handle double counting and PT + EE) for correct integration
        const double weight = mesh.point_area(vi) / 4;
        const Eigen::SparseVector<double> weight_gradient =
            mesh.point_area_gradient(vi) / 4;
#else
        const double weight = 1;
        const Eigen::SparseVector<double> weight_gradient(V.rows());
#endif

        switch (dtype) {
        case PointTriangleDistanceType::P_T0:
            add_vertex_vertex_constraint(
                C_vv, vv_to_index, vi, f0i, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::P_T1:
            add_vertex_vertex_constraint(
                C_vv, vv_to_index, vi, f1i, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::P_T2:
            add_vertex_vertex_constraint(
                C_vv, vv_to_index, vi, f2i, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::P_E0:
            add_edge_vertex_constraint(
                C_ev, ev_to_index, F2E(fi, 0), vi, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::P_E1:
            add_edge_vertex_constraint(
                C_ev, ev_to_index, F2E(fi, 1), vi, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::P_E2:
            add_edge_vertex_constraint(
                C_ev, ev_to_index, F2E(fi, 2), vi, weight, weight_gradient);
            break;

        case PointTriangleDistanceType::P_T:
            C_fv.emplace_back(fi, vi);
            C_fv.back().weight = weight;
#ifdef IPC_TOOLKIT_COMPUTE_SHAPE_DERIVATIVE
            C_fv.back().weight_gradient = weight_gradient;
#endif
            break;
        }
    }
}

void merge_thread_local_constraints(
    const tbb::enumerable_thread_specific<ConstraintBuilder>& local_storage,
    Constraints& constraints)
{
    auto& [C_vv, C_ev, C_ee, C_fv, _] = constraints;
    unordered_map<VertexVertexConstraint, long> vv_to_index;
    unordered_map<EdgeVertexConstraint, long> ev_to_index;

    // size up the hash items
    size_t n_vv = 0, n_ev = 0, n_ee = 0, n_fv = 0;
    for (const auto& storage : local_storage) {
        // This is an conservative estimate
        n_vv += storage.constraint_set.vv_constraints.size();
        n_ev += storage.constraint_set.ev_constraints.size();
        n_ee += storage.constraint_set.ee_constraints.size();
        n_fv += storage.constraint_set.fv_constraints.size();
    }
    C_vv.reserve(n_vv);
    C_ev.reserve(n_ev);
    C_ee.reserve(n_ee);
    C_fv.reserve(n_fv);

    // merge
    for (const auto& storage : local_storage) {
        auto& [lC_vv, lC_ev, lC_ee, lC_fv, __] = storage.constraint_set;

        if (C_vv.empty()) {
            vv_to_index = storage.vv_to_index;
            C_vv.insert(C_vv.end(), lC_vv.begin(), lC_vv.end());
        } else {
            for (const VertexVertexConstraint& vv : lC_vv) {
#ifdef IPC_TOOLKIT_COMPUTE_SHAPE_DERIVATIVE
                const Eigen::SparseVector<double>& weight_gradient =
                    vv.weight_gradient;
#else
                const Eigen::SparseVector<double> weight_gradient;
#endif
                add_vertex_vertex_constraint(
                    C_vv, vv_to_index, vv.vertex0_index, vv.vertex1_index,
                    vv.weight, weight_gradient);
            }
        }

        if (C_ev.empty()) {
            ev_to_index = storage.ev_to_index;
            C_ev.insert(C_ev.end(), lC_ev.begin(), lC_ev.end());
        } else {
            for (const EdgeVertexConstraint& ev : lC_ev) {
#ifdef IPC_TOOLKIT_COMPUTE_SHAPE_DERIVATIVE
                const Eigen::SparseVector<double>& weight_gradient =
                    ev.weight_gradient;
#else
                const Eigen::SparseVector<double> weight_gradient;
#endif
                add_edge_vertex_constraint(
                    C_ev, ev_to_index, ev.edge_index, ev.vertex_index,
                    ev.weight, weight_gradient);
            }
        }

        C_ee.insert(C_ee.end(), lC_ee.begin(), lC_ee.end());
        C_fv.insert(C_fv.end(), lC_fv.begin(), lC_fv.end());
    }
}

} // namespace ipc
