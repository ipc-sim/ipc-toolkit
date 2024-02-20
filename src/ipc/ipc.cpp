#include "ipc.hpp"

#include <ipc/candidates/candidates.hpp>
#include <ipc/utils/intersection.hpp>
#include <ipc/utils/world_bbox_diagonal_length.hpp>

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA
#include <scalable_ccd/cuda/ipc_ccd_strategy.hpp>
#endif

#include <igl/predicates/segment_segment_intersect.h>

namespace ipc {

bool is_step_collision_free(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const BroadPhaseMethod broad_phase_method,
    const double min_distance,
    const double tolerance,
    const long max_iterations)
{
    assert(vertices_t0.rows() == mesh.num_vertices());
    assert(vertices_t1.rows() == mesh.num_vertices());

    // Broad phase
    Candidates candidates;
    candidates.build(
        mesh, vertices_t0, vertices_t1,
        /*inflation_radius=*/min_distance / 2, broad_phase_method);

    // Narrow phase
    return candidates.is_step_collision_free(
        mesh, vertices_t0, vertices_t1, min_distance, tolerance,
        max_iterations);
}

// ============================================================================

double compute_collision_free_stepsize(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const BroadPhaseMethod broad_phase_method,
    const double min_distance,
    const double tolerance,
    const long max_iterations)
{
    assert(vertices_t0.rows() == mesh.num_vertices());
    assert(vertices_t1.rows() == mesh.num_vertices());

    if (broad_phase_method == BroadPhaseMethod::SWEEP_AND_TINIEST_QUEUE) {
#ifdef IPC_TOOLKIT_WITH_CUDA
        if (vertices_t0.cols() != 3) {
            throw std::runtime_error(
                "Sweep and Tiniest Queue is only supported in 3D!");
        }
        // TODO: Use correct min_distance
        const double step_size = scalable_ccd::cuda::ipc_ccd_strategy(
            vertices_t0, vertices_t1, mesh.edges(), mesh.faces(),
            /*min_distance=*/0.0, max_iterations, tolerance);
        if (step_size < 1.0) {
            return 0.8 * step_size;
        }
        return 1.0;
#else
        throw std::runtime_error(
            "GPU Sweep and Tiniest Queue is disabled because CUDA is disabled!");
#endif
    }

    // Broad phase
    Candidates candidates;
    candidates.build(
        mesh, vertices_t0, vertices_t1, /*inflation_radius=*/min_distance / 2,
        broad_phase_method);

    // Narrow phase
    return candidates.compute_collision_free_stepsize(
        mesh, vertices_t0, vertices_t1, min_distance, tolerance,
        max_iterations);
}

// ============================================================================

bool has_intersections(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const BroadPhaseMethod broad_phase_method)
{
    assert(vertices.rows() == mesh.num_vertices());

    const double conservative_inflation_radius =
        1e-6 * world_bbox_diagonal_length(vertices);

    std::shared_ptr<BroadPhase> broad_phase =
        BroadPhase::make_broad_phase(broad_phase_method);
    broad_phase->can_vertices_collide = mesh.can_collide;

    broad_phase->build(
        vertices, mesh.edges(), mesh.faces(), conservative_inflation_radius);

    if (vertices.cols() == 2) {
        // Need to check segment-segment intersections in 2D
        std::vector<EdgeEdgeCandidate> ee_candidates;

        broad_phase->detect_edge_edge_candidates(ee_candidates);
        broad_phase->clear();

        // narrow-phase using igl
        igl::predicates::exactinit();
        for (const auto& [ea_id, eb_id] : ee_candidates) {
            if (igl::predicates::segment_segment_intersect(
                    vertices.row(mesh.edges()(ea_id, 0)).head<2>(),
                    vertices.row(mesh.edges()(ea_id, 1)).head<2>(),
                    vertices.row(mesh.edges()(eb_id, 0)).head<2>(),
                    vertices.row(mesh.edges()(eb_id, 1)).head<2>())) {
                return true;
            }
        }
    } else {
        // Need to check segment-triangle intersections in 3D
        assert(vertices.cols() == 3);

        std::vector<EdgeFaceCandidate> ef_candidates;
        broad_phase->detect_edge_face_candidates(ef_candidates);
        broad_phase->clear();

        for (const auto& [e_id, f_id] : ef_candidates) {
            if (is_edge_intersecting_triangle(
                    vertices.row(mesh.edges()(e_id, 0)),
                    vertices.row(mesh.edges()(e_id, 1)),
                    vertices.row(mesh.faces()(f_id, 0)),
                    vertices.row(mesh.faces()(f_id, 1)),
                    vertices.row(mesh.faces()(f_id, 2)))) {
                return true;
            }
        }
    }

    return false;
}
} // namespace ipc
