#include <ipc/ipc.hpp>

#include <stdexcept> // std::runtime_error
#include <algorithm> // std::min/max

#include <tbb/mutex.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

#include <igl/predicates/segment_segment_intersect.h>

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

/// Find the index of the undirected edge (e0, e1)
long find_edge(const Eigen::MatrixXi& E, long e0, long e1)
{
    for (long i = 0; i < E.rows(); i++) {
        if ((E(i, 0) == e0 && E(i, 1) == e1)
            || (E(i, 1) == e0 && E(i, 0) == e1)) {
            return i;
        }
    }
    throw std::runtime_error("Unable to find edge!");
}

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
    const Eigen::MatrixXd& V_rest,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double dhat,
    Constraints& constraint_set,
    const Eigen::MatrixXi& F2E,
    double dmin,
    const BroadPhaseMethod& method,
    bool ignore_codimensional_vertices,
    const std::function<bool(size_t, size_t)>& can_collide)
{
    double inflation_radius = (dhat + dmin) / 1.99; // Conservative inflation

    Candidates candidates;
    construct_collision_candidates(
        V, E, F, candidates, inflation_radius, method,
        ignore_codimensional_vertices, can_collide);

    construct_constraint_set(
        candidates, V_rest, V, E, F, dhat, constraint_set, F2E, dmin);
}

void construct_constraint_set(
    const Candidates& candidates,
    const Eigen::MatrixXd& V_rest,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double dhat,
    Constraints& constraint_set,
    const Eigen::MatrixXi& F2E,
    double dmin)
{
    constraint_set.clear();

    double dhat_squared = dhat * dhat;
    double dmin_squared = dmin * dmin;

    // Cull the candidates by measuring the distance and dropping those that are
    // greater than dhat.

    // Store the indices to VV and EV pairs to avoid duplicates.
    auto vv_hash = [&V](const VertexVertexConstraint& vv_constraint) -> size_t {
        // The vertex-vertex pair should be order independent
        long min_vi =
            std::min(vv_constraint.vertex0_index, vv_constraint.vertex1_index);
        long max_vi =
            std::max(vv_constraint.vertex0_index, vv_constraint.vertex1_index);
        return size_t(max_vi * V.rows() + min_vi);
    };
    unordered_map<VertexVertexConstraint, long, decltype(vv_hash)> vv_to_index(
        /*min_buckets=*/candidates.size(), vv_hash);

    auto ev_hash = [&V](const EdgeVertexConstraint& ev_constraint) -> size_t {
        // There are max E.rows() * V.rows() constraints
        return size_t(
            ev_constraint.edge_index * V.rows() + ev_constraint.vertex_index);
    };
    unordered_map<EdgeVertexConstraint, long, decltype(ev_hash)> ev_to_index(
        /*min_buckets=*/candidates.size(), ev_hash);

    for (const auto& ev_candidate : candidates.ev_candidates) {
        long vi = ev_candidate.vertex_index;
        long e0i = E(ev_candidate.edge_index, 0);
        long e1i = E(ev_candidate.edge_index, 1);

        PointEdgeDistanceType dtype =
            point_edge_distance_type(V.row(vi), V.row(e0i), V.row(e1i));

        double distance_sqr =
            point_edge_distance(V.row(vi), V.row(e0i), V.row(e1i), dtype);

        if (distance_sqr - dmin_squared < 2 * dmin * dhat + dhat_squared) {
            switch (dtype) {
            case PointEdgeDistanceType::P_E0:
                add_vertex_vertex_constraint(
                    constraint_set.vv_constraints, vv_to_index, vi, e0i);
                break;

            case PointEdgeDistanceType::P_E1:
                add_vertex_vertex_constraint(
                    constraint_set.vv_constraints, vv_to_index, vi, e1i);
                break;

            case PointEdgeDistanceType::P_E:
                // ev_candidates is a set, so no duplicate EV constraints
                constraint_set.ev_constraints.emplace_back(ev_candidate);
                ev_to_index.emplace(
                    constraint_set.ev_constraints.back(),
                    constraint_set.ev_constraints.size() - 1);
                break;
            }
        }
    }

    for (const auto& ee_candidate : candidates.ee_candidates) {
        long eai = ee_candidate.edge0_index, ebi = ee_candidate.edge1_index;
        long ea0i = E(eai, 0), ea1i = E(eai, 1);
        long eb0i = E(ebi, 0), eb1i = E(ebi, 1);

        EdgeEdgeDistanceType dtype = edge_edge_distance_type(
            V.row(ea0i), V.row(ea1i), V.row(eb0i), V.row(eb1i));

        double distance_sqr = edge_edge_distance(
            V.row(ea0i), V.row(ea1i), V.row(eb0i), V.row(eb1i), dtype);

        if (distance_sqr - dmin_squared < 2 * dmin * dhat + dhat_squared) {
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
                    constraint_set.vv_constraints, vv_to_index, ea0i, eb0i);
                break;

            case EdgeEdgeDistanceType::EA0_EB1:
                add_vertex_vertex_constraint(
                    constraint_set.vv_constraints, vv_to_index, ea0i, eb1i);
                break;

            case EdgeEdgeDistanceType::EA1_EB0:
                add_vertex_vertex_constraint(
                    constraint_set.vv_constraints, vv_to_index, ea1i, eb0i);
                break;

            case EdgeEdgeDistanceType::EA1_EB1:
                add_vertex_vertex_constraint(
                    constraint_set.vv_constraints, vv_to_index, ea1i, eb1i);
                break;

            case EdgeEdgeDistanceType::EA_EB0:
                add_edge_vertex_constraint(
                    constraint_set.ev_constraints, ev_to_index, eai, eb0i);
                break;

            case EdgeEdgeDistanceType::EA_EB1:
                add_edge_vertex_constraint(
                    constraint_set.ev_constraints, ev_to_index, eai, eb1i);
                break;

            case EdgeEdgeDistanceType::EA0_EB:
                add_edge_vertex_constraint(
                    constraint_set.ev_constraints, ev_to_index, ebi, ea0i);
                break;

            case EdgeEdgeDistanceType::EA1_EB:
                add_edge_vertex_constraint(
                    constraint_set.ev_constraints, ev_to_index, ebi, ea1i);
                break;

            case EdgeEdgeDistanceType::EA_EB:
                constraint_set.ee_constraints.emplace_back(ee_candidate, eps_x);
                break;
            }
        }
    }

    for (const auto& fv_candidate : candidates.fv_candidates) {
        long vi = fv_candidate.vertex_index;
        long fi = fv_candidate.face_index;
        long f0i = F(fi, 0);
        long f1i = F(fi, 1);
        long f2i = F(fi, 2);

        // Compute distance type
        PointTriangleDistanceType dtype = point_triangle_distance_type(
            V.row(fv_candidate.vertex_index), //
            V.row(f0i), V.row(f1i), V.row(f2i));

        double distance_sqr = point_triangle_distance(
            V.row(fv_candidate.vertex_index), //
            V.row(f0i), V.row(f1i), V.row(f2i), dtype);

        if (distance_sqr - dmin_squared < 2 * dmin * dhat + dhat_squared) {
            switch (dtype) {
            case PointTriangleDistanceType::P_T0:
                add_vertex_vertex_constraint(
                    constraint_set.vv_constraints, vv_to_index, vi, f0i);
                break;

            case PointTriangleDistanceType::P_T1:
                add_vertex_vertex_constraint(
                    constraint_set.vv_constraints, vv_to_index, vi, f1i);
                break;

            case PointTriangleDistanceType::P_T2:
                add_vertex_vertex_constraint(
                    constraint_set.vv_constraints, vv_to_index, vi, f2i);
                break;

            case PointTriangleDistanceType::P_E0:
                add_edge_vertex_constraint(
                    constraint_set.ev_constraints, ev_to_index,
                    F2E.rows() > fi ? F2E(fi, 0) : find_edge(E, f0i, f1i), vi);
                break;

            case PointTriangleDistanceType::P_E1:
                add_edge_vertex_constraint(
                    constraint_set.ev_constraints, ev_to_index,
                    F2E.rows() > fi ? F2E(fi, 1) : find_edge(E, f1i, f2i), vi);
                break;

            case PointTriangleDistanceType::P_E2:
                add_edge_vertex_constraint(
                    constraint_set.ev_constraints, ev_to_index,
                    F2E.rows() > fi ? F2E(fi, 2) : find_edge(E, f2i, f0i), vi);
                break;

            case PointTriangleDistanceType::P_T:
                constraint_set.fv_constraints.emplace_back(fv_candidate);
                break;
            }
        }
    }

    for (size_t ci = 0; ci < constraint_set.size(); ci++) {
        constraint_set[ci].minimum_distance = dmin;
    }
}

///////////////////////////////////////////////////////////////////////////////

double compute_barrier_potential(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Constraints& constraint_set,
    double dhat)
{
    if (constraint_set.empty()) {
        return 0;
    }

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
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Constraints& constraint_set,
    double dhat)
{
    if (constraint_set.empty()) {
        return Eigen::VectorXd::Zero(V.size());
    }

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
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Constraints& constraint_set,
    double dhat,
    bool project_hessian_to_psd)
{
    if (constraint_set.empty()) {
        return Eigen::SparseMatrix<double>(V.size(), V.size());
    }

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
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const BroadPhaseMethod& method,
    double tolerance,
    int max_iterations,
    bool ignore_codimensional_vertices,
    const std::function<bool(size_t, size_t)>& can_collide)
{
    // Broad phase
    Candidates candidates;
    construct_collision_candidates(
        V0, V1, E, F, candidates, /*inflation_radius=*/0, method,
        ignore_codimensional_vertices, can_collide);

    // Narrow phase
    return is_step_collision_free(
        candidates, V0, V1, E, F, tolerance, max_iterations);
}

bool is_step_collision_free(
    const Candidates& candidates,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double tolerance,
    int max_iterations)
{
    assert(V0.cols() == V1.cols());

    // Narrow phase
    for (const auto& ev_candidate : candidates.ev_candidates) {
        double toi;
        bool is_collision = point_edge_ccd(
            // Point at t=0
            V0.row(ev_candidate.vertex_index),
            // Edge at t=0
            V0.row(E(ev_candidate.edge_index, 0)),
            V0.row(E(ev_candidate.edge_index, 1)),
            // Point at t=1
            V1.row(ev_candidate.vertex_index),
            // Edge at t=1
            V1.row(E(ev_candidate.edge_index, 0)),
            V1.row(E(ev_candidate.edge_index, 1)), //
            toi, /*tmax=*/1.0, tolerance, max_iterations);

        if (is_collision) {
            return false;
        }
    }

    for (const auto& ee_candidate : candidates.ee_candidates) {
        double toi;
        bool is_collision = edge_edge_ccd(
            // Edge 1 at t=0
            V0.row(E(ee_candidate.edge0_index, 0)),
            V0.row(E(ee_candidate.edge0_index, 1)),
            // Edge 2 at t=0
            V0.row(E(ee_candidate.edge1_index, 0)),
            V0.row(E(ee_candidate.edge1_index, 1)),
            // Edge 1 at t=1
            V1.row(E(ee_candidate.edge0_index, 0)),
            V1.row(E(ee_candidate.edge0_index, 1)),
            // Edge 2 at t=1
            V1.row(E(ee_candidate.edge1_index, 0)),
            V1.row(E(ee_candidate.edge1_index, 1)), //
            toi, /*tmax=*/1.0, tolerance, max_iterations);

        if (is_collision) {
            return false;
        }
    }

    for (const auto& fv_candidate : candidates.fv_candidates) {
        double toi;
        bool is_collision = point_triangle_ccd(
            // Point at t=0
            V0.row(fv_candidate.vertex_index),
            // Triangle at t = 0
            V0.row(F(fv_candidate.face_index, 0)),
            V0.row(F(fv_candidate.face_index, 1)),
            V0.row(F(fv_candidate.face_index, 2)),
            // Point at t=1
            V1.row(fv_candidate.vertex_index),
            // Triangle at t = 1
            V1.row(F(fv_candidate.face_index, 0)),
            V1.row(F(fv_candidate.face_index, 1)),
            V1.row(F(fv_candidate.face_index, 2)), //
            toi, /*tmax=*/1.0, tolerance, max_iterations);

        if (is_collision) {
            return false;
        }
    }

    return true;
}

///////////////////////////////////////////////////////////////////////////////

double compute_collision_free_stepsize(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const BroadPhaseMethod& method,
    double tolerance,
    int max_iterations,
    bool ignore_codimensional_vertices,
    const std::function<bool(size_t, size_t)>& can_collide)
{
    // Broad phase
    Candidates candidates;
    construct_collision_candidates(
        V0, V1, E, F, candidates, /*inflation_radius=*/0, method,
        ignore_codimensional_vertices, can_collide);

    // Narrow phase
    return compute_collision_free_stepsize(
        candidates, V0, V1, E, F, tolerance, max_iterations);
}

double compute_collision_free_stepsize(
    const Candidates& candidates,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double tolerance,
    int max_iterations)
{
    assert(V0.cols() == V1.cols());

    if (candidates.empty()) {
        return 1; // No possible collisions, so can take full step.
    }

    // Narrow phase
    double earliest_toi = 1;
    tbb::mutex earliest_toi_mutex;

    const size_t num_ev = candidates.ev_candidates.size();
    const size_t num_ee = candidates.ee_candidates.size();
    const size_t num_fv = candidates.fv_candidates.size();

    // Do a single block range over all three candidate vectors
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, candidates.size()),
        [&](tbb::blocked_range<size_t> r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                double toi = std::numeric_limits<double>::infinity();
                bool are_colliding;

                size_t ci = i;
                if (ci < num_ev) {
                    const auto& ev_candidate = candidates.ev_candidates[ci];
                    const auto& ei = ev_candidate.edge_index;
                    const auto& vi = ev_candidate.vertex_index;
                    are_colliding = point_edge_ccd(
                        // Point at t=0
                        V0.row(vi),
                        // Edge at t=0
                        V0.row(E(ei, 0)), V0.row(E(ei, 1)),
                        // Point at t=1
                        V1.row(vi),
                        // Edge at t=1
                        V1.row(E(ei, 0)), V1.row(E(ei, 1)), //
                        toi, /*tmax=*/earliest_toi, tolerance, max_iterations);
                } else if ((ci -= num_ev) < num_ee) {
                    const auto& ee_candidate = candidates.ee_candidates[ci];
                    const auto& eai = ee_candidate.edge0_index;
                    const auto& ebi = ee_candidate.edge1_index;
                    are_colliding = edge_edge_ccd(
                        // Edge 1 at t=0
                        V0.row(E(eai, 0)), V0.row(E(eai, 1)),
                        // Edge 2 at t=0
                        V0.row(E(ebi, 0)), V0.row(E(ebi, 1)),
                        // Edge 1 at t=1
                        V1.row(E(eai, 0)), V1.row(E(eai, 1)),
                        // Edge 2 at t=1
                        V1.row(E(ebi, 0)), V1.row(E(ebi, 1)), //
                        toi, /*tmax=*/earliest_toi, tolerance, max_iterations);
                } else {
                    ci -= num_ee;
                    assert(ci < num_fv);
                    const auto& fv_candidate = candidates.fv_candidates[ci];
                    const auto& fi = fv_candidate.face_index;
                    const auto& vi = fv_candidate.vertex_index;
                    are_colliding = point_triangle_ccd(
                        // Point at t=0
                        V0.row(vi),
                        // Triangle at t = 0
                        V0.row(F(fi, 0)), V0.row(F(fi, 1)), V0.row(F(fi, 2)),
                        // Point at t=1
                        V1.row(vi),
                        // Triangle at t = 1
                        V1.row(F(fi, 0)), V1.row(F(fi, 1)), V1.row(F(fi, 2)),
                        toi, /*tmax=*/earliest_toi, tolerance, max_iterations);
                }

                if (are_colliding) {
                    tbb::mutex::scoped_lock lock(earliest_toi_mutex);
                    if (toi < earliest_toi) {
                        earliest_toi = toi;
                    }
                }
            }
        });

    assert(earliest_toi >= 0 && earliest_toi <= 1.0);
    return earliest_toi;
}

///////////////////////////////////////////////////////////////////////////////

// NOTE: Actually distance squared
double compute_minimum_distance(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Constraints& constraint_set)
{
    if (constraint_set.empty()) {
        return std::numeric_limits<double>::infinity();
    }

    tbb::enumerable_thread_specific<double> storage(
        std::numeric_limits<double>::infinity());

    const size_t num_vv = constraint_set.vv_constraints.size();
    const size_t num_ev = constraint_set.ev_constraints.size();
    const size_t num_ee = constraint_set.ee_constraints.size();
    const size_t num_fv = constraint_set.fv_constraints.size();

    // Do a single block range over all constraint vectors
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, constraint_set.size()),
        [&](tbb::blocked_range<size_t> r) {
            auto& local_min_dist = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {

                double dist;
                size_t ci = i;
                if (ci < num_vv) {
                    const auto& vv_constraint =
                        constraint_set.vv_constraints[ci];
                    const auto& vai = vv_constraint.vertex0_index;
                    const auto& vbi = vv_constraint.vertex1_index;
                    dist = point_point_distance(V.row(vai), V.row(vbi));
                } else if ((ci -= num_vv) < num_ev) {
                    const auto& ev_constraint =
                        constraint_set.ev_constraints[ci];
                    const auto& ei = ev_constraint.edge_index;
                    const auto& vi = ev_constraint.vertex_index;
                    dist = point_edge_distance(
                        V.row(vi), V.row(E(ei, 0)), V.row(E(ei, 1)),
                        PointEdgeDistanceType::P_E);
                } else if ((ci -= num_ev) < num_ee) {
                    const auto& ee_constraint =
                        constraint_set.ee_constraints[ci];
                    const auto& eai = ee_constraint.edge0_index;
                    const auto& ebi = ee_constraint.edge1_index;
                    // The distance type is unknown because of mollified PP and
                    // PE constraints where also added as EE constraints.
                    dist = edge_edge_distance(
                        V.row(E(eai, 0)), V.row(E(eai, 1)), //
                        V.row(E(ebi, 0)), V.row(E(ebi, 1)));
                } else {
                    ci -= num_ee;
                    assert(ci < num_fv);
                    const auto& fv_constraint =
                        constraint_set.fv_constraints[ci];
                    const auto& fi = fv_constraint.face_index;
                    const auto& vi = fv_constraint.vertex_index;
                    dist = point_triangle_distance(
                        V.row(vi), V.row(F(fi, 0)), V.row(F(fi, 1)),
                        V.row(F(fi, 2)), PointTriangleDistanceType::P_T);
                }

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
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const std::function<bool(size_t, size_t)>& can_collide)
{
    // TODO: Expose the broad-phase method
    HashGrid hash_grid;
    double conservative_inflation_radius = 1e-2 * world_bbox_diagonal_length(V);
    hash_grid.resize(V, E, conservative_inflation_radius);
    hash_grid.addEdges(V, E, conservative_inflation_radius);
    if (V.cols() == 3) {
        // These are not needed for 2D
        hash_grid.addFaces(V, F, conservative_inflation_radius);
    }

    if (V.cols() == 2) { // Need to check segment-segment intersections in 2D
        std::vector<EdgeEdgeCandidate> ee_candidates;
        hash_grid.getEdgeEdgePairs(E, ee_candidates, can_collide);

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
        hash_grid.getEdgeFacePairs(E, F, ef_candidates, can_collide);

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
