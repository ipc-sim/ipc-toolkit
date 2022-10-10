#include "ipc.hpp"

#include <ipc/ccd/ccd.hpp>
#include <ipc/broad_phase/hash_grid.hpp>
#include <ipc/utils/local_to_global.hpp>
#include <ipc/utils/intersection.hpp>
#include <ipc/utils/world_bbox_diagonal_length.hpp>

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA
#include <ccdgpu/helper.cuh>
#endif

#include <igl/predicates/segment_segment_intersect.h>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

#define IPC_EARLIEST_TOI_USE_MUTEX
#ifdef IPC_EARLIEST_TOI_USE_MUTEX
#include <mutex>
#endif

#include <algorithm> // std::min/max

namespace ipc {

bool is_step_collision_free(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const BroadPhaseMethod method,
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
    const BroadPhaseMethod method,
    const double tolerance,
    const long max_iterations)
{
    assert(V0.rows() == mesh.num_vertices());
    assert(V1.rows() == mesh.num_vertices());
    const Eigen::MatrixXi& E = mesh.edges();
    const Eigen::MatrixXi& F = mesh.faces();

    if (method == BroadPhaseMethod::SWEEP_AND_TINIEST_QUEUE_GPU) {
#ifdef IPC_TOOLKIT_WITH_CUDA
        double min_distance = 0; // TODO
        const double step_size = ccd::gpu::compute_toi_strategy(
            V0, V1, E, F, max_iterations, min_distance, tolerance);
        if (step_size < 1.0) {
            return 0.8 * step_size;
        }
        return 1.0;
#else
        throw std::runtime_error("GPU Sweep and Tiniest Queue is disabled "
                                 "because CUDA is disabled!");
#endif
    }

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

bool has_intersections(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const BroadPhaseMethod method)
{
    assert(V.rows() == mesh.num_vertices());
    const Eigen::MatrixXi& E = mesh.edges();
    const Eigen::MatrixXi& F = mesh.faces();

    const double conservative_inflation_radius =
        1e-6 * world_bbox_diagonal_length(V);

    // TODO: Expose the broad-phase method
    std::unique_ptr<BroadPhase> broad_phase =
        BroadPhase::make_broad_phase(method);
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
