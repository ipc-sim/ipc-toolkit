#include <ipc/ipc.hpp>

#include <stdexcept> // std::runtime_error
#include <algorithm> // std::min/max

#define IPC_EARLIEST_TOI_USE_MUTEX
#ifdef IPC_EARLIEST_TOI_USE_MUTEX
#include <mutex>
#endif
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

#include <igl/predicates/segment_segment_intersect.h>
#ifdef IPC_TOOLKIT_WITH_CUDA
#include <ccdgpu/helper.cuh>
#endif

#include <ipc/ccd/ccd.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/broad_phase/hash_grid.hpp>
#include <ipc/utils/local_to_global.hpp>
#include <ipc/utils/intersection.hpp>
#include <ipc/utils/world_bbox_diagonal_length.hpp>
#include <ipc/utils/unordered_map_and_set.hpp>

namespace ipc {

template <typename Hash>
void add_vertex_vertex_constraint(
    std::vector<VertexVertexConstraint>& vv_constraints,
    unordered_map<VertexVertexConstraint, long, Hash>& vv_to_index,
    const long v0i,
    const long v1i)
{
    VertexVertexConstraint vv_constraint(v0i, v1i);
    auto found_item = vv_to_index.find(vv_constraint);
    if (found_item != vv_to_index.end()) {
        // Constraint already exists, so increase multiplicity
        vv_constraints[found_item->second].multiplicity++;
    } else {
        // New constraint, so add it to the end of vv_constraints
        vv_to_index.emplace(vv_constraint, vv_constraints.size());
        vv_constraints.push_back(vv_constraint);
    }
}

template <typename Hash>
void add_edge_vertex_constraint(
    std::vector<EdgeVertexConstraint>& ev_constraints,
    unordered_map<EdgeVertexConstraint, long, Hash>& ev_to_index,
    const long ei,
    const long vi)
{
    EdgeVertexConstraint ev_constraint(ei, vi);
    auto found_item = ev_to_index.find(ev_constraint);
    if (found_item != ev_to_index.end()) {
        // Constraint already exists, so increase multiplicity
        ev_constraints[found_item->second].multiplicity++;
    } else {
        // New constraint, so add it to the end of vv_constraints
        ev_to_index.emplace(ev_constraint, ev_constraints.size());
        ev_constraints.push_back(ev_constraint);
    }
}

void construct_constraint_set(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const double dhat,
    Constraints& constraint_set,
    const double dmin,
    const BroadPhaseMethod& method)
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

    const Eigen::MatrixXd& V_rest = mesh.vertices_at_rest();
    const Eigen::MatrixXi& E = mesh.edges();
    const Eigen::MatrixXi& F = mesh.faces();
    const Eigen::MatrixXi& F2E = mesh.faces_to_edges();

    // Cull the candidates by measuring the distance and dropping those that are
    // greater than dhat.
    const double offset_sqr = (dmin + dhat) * (dmin + dhat);
    auto is_active = [&](double distance_sqr) {
        return distance_sqr < offset_sqr;
    };

    // Store the indices to VV and EV pairs to avoid duplicates.
    unordered_map<VertexVertexConstraint, long> vv_to_index;
    unordered_map<EdgeVertexConstraint, long> ev_to_index;

    std::mutex vv_mutex, ev_mutex, ee_mutex, fv_mutex;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), candidates.ev_candidates.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                const auto& ev_candidate = candidates.ev_candidates[i];
                const auto& [ei, vi] = ev_candidate;
                long e0i = E(ei, 0), e1i = E(ei, 1);

                PointEdgeDistanceType dtype =
                    point_edge_distance_type(V.row(vi), V.row(e0i), V.row(e1i));

                double distance_sqr = point_edge_distance(
                    V.row(vi), V.row(e0i), V.row(e1i), dtype);

                if (is_active(distance_sqr)) {
                    switch (dtype) {
                    case PointEdgeDistanceType::P_E0: {
                        std::lock_guard<std::mutex> lock(vv_mutex);
                        add_vertex_vertex_constraint(
                            constraint_set.vv_constraints, vv_to_index, vi,
                            e0i);
                    } break;

                    case PointEdgeDistanceType::P_E1: {
                        std::lock_guard<std::mutex> lock(vv_mutex);
                        add_vertex_vertex_constraint(
                            constraint_set.vv_constraints, vv_to_index, vi,
                            e1i);
                    } break;

                    case PointEdgeDistanceType::P_E: {
                        std::lock_guard<std::mutex> lock(ev_mutex);
                        // ev_candidates is a set, so no duplicate EV
                        // constraints
                        constraint_set.ev_constraints.emplace_back(
                            ev_candidate);
                        ev_to_index.emplace(
                            constraint_set.ev_constraints.back(),
                            constraint_set.ev_constraints.size() - 1);
                    } break;
                    }
                }
            }
        });

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), candidates.ee_candidates.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                const auto& ee_candidate = candidates.ee_candidates[i];
                const auto& [eai, ebi] = ee_candidate;
                long ea0i = E(eai, 0), ea1i = E(eai, 1);
                long eb0i = E(ebi, 0), eb1i = E(ebi, 1);

                EdgeEdgeDistanceType dtype = edge_edge_distance_type(
                    V.row(ea0i), V.row(ea1i), V.row(eb0i), V.row(eb1i));

                double distance_sqr = edge_edge_distance(
                    V.row(ea0i), V.row(ea1i), V.row(eb0i), V.row(eb1i), dtype);

                if (is_active(distance_sqr)) {
                    double eps_x = edge_edge_mollifier_threshold(
                        V_rest.row(ea0i), V_rest.row(ea1i), //
                        V_rest.row(eb0i), V_rest.row(eb1i));
                    double ee_cross_norm_sqr = edge_edge_cross_squarednorm(
                        V.row(ea0i), V.row(ea1i), V.row(eb0i), V.row(eb1i));
                    if (ee_cross_norm_sqr < eps_x) {
                        // NOTE: This may not actually be the distance type, but
                        // all EE pairs requiring mollification must be
                        // mollified later.
                        dtype = EdgeEdgeDistanceType::EA_EB;
                    }

                    switch (dtype) {
                    case EdgeEdgeDistanceType::EA0_EB0: {
                        std::lock_guard<std::mutex> lock(vv_mutex);
                        add_vertex_vertex_constraint(
                            constraint_set.vv_constraints, vv_to_index, ea0i,
                            eb0i);
                    } break;

                    case EdgeEdgeDistanceType::EA0_EB1: {
                        std::lock_guard<std::mutex> lock(vv_mutex);
                        add_vertex_vertex_constraint(
                            constraint_set.vv_constraints, vv_to_index, ea0i,
                            eb1i);
                    } break;

                    case EdgeEdgeDistanceType::EA1_EB0: {
                        std::lock_guard<std::mutex> lock(vv_mutex);
                        add_vertex_vertex_constraint(
                            constraint_set.vv_constraints, vv_to_index, ea1i,
                            eb0i);
                    } break;

                    case EdgeEdgeDistanceType::EA1_EB1: {
                        std::lock_guard<std::mutex> lock(vv_mutex);
                        add_vertex_vertex_constraint(
                            constraint_set.vv_constraints, vv_to_index, ea1i,
                            eb1i);
                    } break;

                    case EdgeEdgeDistanceType::EA_EB0: {
                        std::lock_guard<std::mutex> lock(ev_mutex);
                        add_edge_vertex_constraint(
                            constraint_set.ev_constraints, ev_to_index, eai,
                            eb0i);
                    } break;

                    case EdgeEdgeDistanceType::EA_EB1: {
                        std::lock_guard<std::mutex> lock(ev_mutex);
                        add_edge_vertex_constraint(
                            constraint_set.ev_constraints, ev_to_index, eai,
                            eb1i);
                    } break;

                    case EdgeEdgeDistanceType::EA0_EB: {
                        std::lock_guard<std::mutex> lock(ev_mutex);
                        add_edge_vertex_constraint(
                            constraint_set.ev_constraints, ev_to_index, ebi,
                            ea0i);
                    } break;

                    case EdgeEdgeDistanceType::EA1_EB: {
                        std::lock_guard<std::mutex> lock(ev_mutex);
                        add_edge_vertex_constraint(
                            constraint_set.ev_constraints, ev_to_index, ebi,
                            ea1i);
                    } break;

                    case EdgeEdgeDistanceType::EA_EB: {
                        std::lock_guard<std::mutex> lock(ee_mutex);
                        constraint_set.ee_constraints.emplace_back(
                            ee_candidate, eps_x);
                    } break;
                    }
                }
            }
        });

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), candidates.fv_candidates.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                const auto& fv_candidate = candidates.fv_candidates[i];
                const auto& [fi, vi] = fv_candidate;
                long f0i = F(fi, 0), f1i = F(fi, 1), f2i = F(fi, 2);

                // Compute distance type
                PointTriangleDistanceType dtype = point_triangle_distance_type(
                    V.row(fv_candidate.vertex_index), //
                    V.row(f0i), V.row(f1i), V.row(f2i));

                double distance_sqr = point_triangle_distance(
                    V.row(fv_candidate.vertex_index), //
                    V.row(f0i), V.row(f1i), V.row(f2i), dtype);

                if (is_active(distance_sqr)) {
                    switch (dtype) {
                    case PointTriangleDistanceType::P_T0: {
                        std::lock_guard<std::mutex> lock(vv_mutex);
                        add_vertex_vertex_constraint(
                            constraint_set.vv_constraints, vv_to_index, vi,
                            f0i);
                    } break;

                    case PointTriangleDistanceType::P_T1: {
                        std::lock_guard<std::mutex> lock(vv_mutex);
                        add_vertex_vertex_constraint(
                            constraint_set.vv_constraints, vv_to_index, vi,
                            f1i);
                    } break;

                    case PointTriangleDistanceType::P_T2: {
                        std::lock_guard<std::mutex> lock(vv_mutex);
                        add_vertex_vertex_constraint(
                            constraint_set.vv_constraints, vv_to_index, vi,
                            f2i);
                    } break;

                    case PointTriangleDistanceType::P_E0: {
                        std::lock_guard<std::mutex> lock(ev_mutex);
                        add_edge_vertex_constraint(
                            constraint_set.ev_constraints, ev_to_index,
                            F2E(fi, 0), vi);
                    } break;

                    case PointTriangleDistanceType::P_E1: {
                        std::lock_guard<std::mutex> lock(ev_mutex);
                        add_edge_vertex_constraint(
                            constraint_set.ev_constraints, ev_to_index,
                            F2E(fi, 1), vi);
                    } break;

                    case PointTriangleDistanceType::P_E2: {
                        std::lock_guard<std::mutex> lock(ev_mutex);
                        add_edge_vertex_constraint(
                            constraint_set.ev_constraints, ev_to_index,
                            F2E(fi, 2), vi);
                    } break;

                    case PointTriangleDistanceType::P_T: {
                        std::lock_guard<std::mutex> lock(fv_mutex);
                        constraint_set.fv_constraints.emplace_back(
                            fv_candidate);
                    } break;
                    }
                }
            }
        });

    for (size_t ci = 0; ci < constraint_set.size(); ci++) {
        constraint_set[ci].minimum_distance = dmin;
    }
}

///////////////////////////////////////////////////////////////////////////////

double compute_barrier_potential(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const Constraints& constraint_set,
    const double dhat)
{
    assert(V.rows() == mesh.num_vertices());

    if (constraint_set.empty()) {
        return 0;
    }

    const Eigen::MatrixXi& E = mesh.edges();
    const Eigen::MatrixXi& F = mesh.faces();

    tbb::enumerable_thread_specific<double> storage(0);

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), constraint_set.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& local_potential = storage.local();
            for (size_t i = r.begin(); i < r.end(); i++) {
                local_potential +=
                    constraint_set[i].compute_potential(V, E, F, dhat);
            }
        });

    double potential = 0;
    for (const auto& local_potential : storage) {
        potential += local_potential;
    }
    return potential;
}

Eigen::VectorXd compute_barrier_potential_gradient(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const Constraints& constraint_set,
    const double dhat)
{
    assert(V.rows() == mesh.num_vertices());

    if (constraint_set.empty()) {
        return Eigen::VectorXd::Zero(V.size());
    }

    const Eigen::MatrixXi& E = mesh.edges();
    const Eigen::MatrixXi& F = mesh.faces();

    int dim = V.cols();

    tbb::enumerable_thread_specific<Eigen::VectorXd> storage(
        Eigen::VectorXd::Zero(V.size()));

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), constraint_set.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& local_grad = storage.local();
            for (size_t i = r.begin(); i < r.end(); i++) {
                local_gradient_to_global_gradient(
                    constraint_set[i].compute_potential_gradient(V, E, F, dhat),
                    constraint_set[i].vertex_indices(E, F), dim, local_grad);
            }
        });

    Eigen::VectorXd grad = Eigen::VectorXd::Zero(V.size());
    for (const auto& local_grad : storage) {
        grad += local_grad;
    }
    return grad;
}

Eigen::SparseMatrix<double> compute_barrier_potential_hessian(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const Constraints& constraint_set,
    const double dhat,
    const bool project_hessian_to_psd)
{
    assert(V.rows() == mesh.num_vertices());

    if (constraint_set.empty()) {
        return Eigen::SparseMatrix<double>(V.size(), V.size());
    }

    const Eigen::MatrixXi& E = mesh.edges();
    const Eigen::MatrixXi& F = mesh.faces();

    int dim = V.cols();

    tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<double>>>
        storage;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), constraint_set.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& local_hess_triplets = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                local_hessian_to_global_triplets(
                    constraint_set[i].compute_potential_hessian(
                        V, E, F, dhat, project_hessian_to_psd),
                    constraint_set[i].vertex_indices(E, F), dim,
                    local_hess_triplets);
            }
        });

    Eigen::SparseMatrix<double> hess(V.size(), V.size());
    for (const auto& local_hess_triplets : storage) {
        Eigen::SparseMatrix<double> local_hess(V.size(), V.size());
        local_hess.setFromTriplets(
            local_hess_triplets.begin(), local_hess_triplets.end());
        hess += local_hess;
    }
    return hess;
}

///////////////////////////////////////////////////////////////////////////////

bool is_step_collision_free(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const BroadPhaseMethod& method,
    const double tolerance,
    const long max_iterations)
{
    assert(V0.rows() == mesh.num_vertices());
    assert(V1.rows() == mesh.num_vertices());

    // Broad phase
    Candidates candidates;
    construct_collision_candidates(
        mesh, V0, V1, candidates, /*inflation_radius=*/0, method);

    // Narrow phase
    return is_step_collision_free(
        candidates, mesh, V0, V1, tolerance, max_iterations);
}

bool is_step_collision_free(
    const Candidates& candidates,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const double tolerance,
    const long max_iterations)
{
    assert(V0.rows() == mesh.num_vertices());
    assert(V1.rows() == mesh.num_vertices());

    const Eigen::MatrixXi& E = mesh.edges();
    const Eigen::MatrixXi& F = mesh.faces();

    // Narrow phase
    for (size_t i = 0; i < candidates.size(); i++) {
        double toi;
        bool is_collision = candidates[i].ccd(
            V0, V1, E, F, toi, /*tmax=*/1.0, tolerance, max_iterations);

        if (is_collision) {
            return false;
        }
    }

    return true;
}

///////////////////////////////////////////////////////////////////////////////

double compute_collision_free_stepsize(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const BroadPhaseMethod& method,
    const double tolerance,
    const long max_iterations)
{
    assert(V0.rows() == mesh.num_vertices());
    assert(V1.rows() == mesh.num_vertices());
    const Eigen::MatrixXi& E = mesh.edges();
    const Eigen::MatrixXi& F = mesh.faces();

#ifdef IPC_TOOLKIT_WITH_CUDA
    if (method == BroadPhaseMethod::SWEEP_AND_TINIEST_QUEUE_GPU) {
        double min_distance = 0; // TODO
        const double step_size = ccd::gpu::compute_toi_strategy(
            V0, V1, E, F, max_iterations, min_distance, tolerance);
        if (step_size < 1.0) {
            return 0.8 * step_size;
        }
        return 1.0;
    }
#endif

    // Broad phase
    Candidates candidates;
    construct_collision_candidates(
        mesh, V0, V1, candidates, /*inflation_radius=*/0, method);

    // Narrow phase
    double step_size = compute_collision_free_stepsize(
        candidates, mesh, V0, V1, tolerance, max_iterations);

    return step_size;
}

double compute_collision_free_stepsize(
    const Candidates& candidates,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const double tolerance,
    const long max_iterations)
{
    assert(V0.rows() == mesh.num_vertices());
    assert(V1.rows() == mesh.num_vertices());
    const Eigen::MatrixXi& E = mesh.edges();
    const Eigen::MatrixXi& F = mesh.faces();

    if (candidates.empty()) {
        return 1; // No possible collisions, so can take full step.
    }

    // Narrow phase
#ifdef IPC_EARLIEST_TOI_USE_MUTEX
    double earliest_toi = 1;
    std::mutex earliest_toi_mutex;
#else
    tbb::enumerable_thread_specific<double> storage(1);
#endif

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, candidates.size()),
        [&](tbb::blocked_range<size_t> r) {
#ifndef IPC_EARLIEST_TOI_USE_MUTEX
            double& earliest_toi = storage.local();
#endif
            for (size_t i = r.begin(); i < r.end(); i++) {
                double toi = std::numeric_limits<double>::infinity();
                bool are_colliding = candidates[i].ccd(
                    V0, V1, E, F, toi, /*tmax=*/earliest_toi, tolerance,
                    max_iterations);

                if (are_colliding) {
#ifdef IPC_EARLIEST_TOI_USE_MUTEX
                    std::lock_guard<std::mutex> lock(earliest_toi_mutex);
#endif
                    if (toi < earliest_toi) {
                        earliest_toi = toi;
                    }
                }
            }
        });

#ifndef IPC_EARLIEST_TOI_USE_MUTEX
    double earliest_toi = 1;
    for (const auto& local_earliest_toi : storage) {
        earliest_toi = std::min(earliest_toi, local_earliest_toi);
    }
#endif
    assert(earliest_toi >= 0 && earliest_toi <= 1.0);
    return earliest_toi;
}

///////////////////////////////////////////////////////////////////////////////

// NOTE: Actually distance squared
double compute_minimum_distance(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const Constraints& constraint_set)
{
    assert(V.rows() == mesh.num_vertices());

    if (constraint_set.empty()) {
        return std::numeric_limits<double>::infinity();
    }

    const Eigen::MatrixXi& E = mesh.edges();
    const Eigen::MatrixXi& F = mesh.faces();

    tbb::enumerable_thread_specific<double> storage(
        std::numeric_limits<double>::infinity());

    // Do a single block range over all constraint vectors
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, constraint_set.size()),
        [&](tbb::blocked_range<size_t> r) {
            double& local_min_dist = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                const double dist = constraint_set[i].compute_distance(V, E, F);

                if (dist < local_min_dist) {
                    local_min_dist = dist;
                }
            }
        });

    double min_dist = std::numeric_limits<double>::infinity();
    for (const auto& local_min_dist : storage) {
        min_dist = std::min(min_dist, local_min_dist);
    }
    return min_dist;
}

///////////////////////////////////////////////////////////////////////////////

bool has_intersections(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const BroadPhaseMethod& method)
{
    assert(V.rows() == mesh.num_vertices());
    const Eigen::MatrixXi& E = mesh.edges();
    const Eigen::MatrixXi& F = mesh.faces();

    const double conservative_inflation_radius = 1e-6
        * world_bbox_diagonal_length(V)

        // TODO: Expose the broad-phase method
        std::unique_ptr<BroadPhase>
            broad_phase = BroadPhase::make_broad_phase(method);
    broad_phase->can_vertices_collide = mesh.can_collide;

    broad_phase->build(V, E, F, conservative_inflation_radius);

    if (V.cols() == 2) { // Need to check segment-segment intersections in 2D
        std::vector<EdgeEdgeCandidate> ee_candidates;

        broad_phase->detect_edge_edge_candidates(ee_candidates);
        broad_phase->clear();

        // narrow-phase using igl
        igl::predicates::exactinit();
        for (const EdgeEdgeCandidate& ee_candidate : ee_candidates) {
            if (igl::predicates::segment_segment_intersect(
                    V.row(E(ee_candidate.edge0_index, 0)).head<2>(),
                    V.row(E(ee_candidate.edge0_index, 1)).head<2>(),
                    V.row(E(ee_candidate.edge1_index, 0)).head<2>(),
                    V.row(E(ee_candidate.edge1_index, 1)).head<2>())) {
                return true;
            }
        }
    } else { // Need to check segment-triangle intersections in 3D
        assert(V.cols() == 3);

        std::vector<EdgeFaceCandidate> ef_candidates;
        broad_phase->detect_edge_face_candidates(ef_candidates);
        broad_phase->clear();

        for (const EdgeFaceCandidate& ef_candidate : ef_candidates) {
            if (is_edge_intersecting_triangle(
                    V.row(E(ef_candidate.edge_index, 0)),
                    V.row(E(ef_candidate.edge_index, 1)),
                    V.row(F(ef_candidate.face_index, 0)),
                    V.row(F(ef_candidate.face_index, 1)),
                    V.row(F(ef_candidate.face_index, 2)))) {
                return true;
            }
        }
    }

    return false;
}
} // namespace ipc
