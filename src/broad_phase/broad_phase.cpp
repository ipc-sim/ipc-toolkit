#include <ipc/broad_phase/broad_phase.hpp>

#include <ipc/broad_phase/brute_force.hpp>
#include <ipc/broad_phase/spatial_hash.hpp>
#include <ipc/broad_phase/hash_grid.hpp>

namespace ipc {

void construct_collision_candidates(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    Candidates& candidates,
    double inflation_radius,
    const BroadPhaseMethod& method,
    bool ignore_codimensional_vertices,
    const std::function<bool(size_t, size_t)>& can_collide)
{
    int dim = V.cols();

    candidates.clear();

    switch (method) {
    case BroadPhaseMethod::BRUTE_FORCE:
        // TODO: Use ignore_codimensional_vertices
        detect_collision_candidates_brute_force(
            V, E, F, candidates,
            /*detect_edge_vertex=*/dim == 2,
            /*detect_edge_edge=*/dim == 3,
            /*detect_face_vertex=*/dim == 3,
            /*perform_aabb_check=*/true,
            /*aabb_inflation_radius=*/inflation_radius,
            ignore_codimensional_vertices, can_collide);
        break;
    case BroadPhaseMethod::HASH_GRID: {
        HashGrid hash_grid;
        hash_grid.resize(V, E, inflation_radius);

        // Assumes the edges connect to all boundary vertices
        if (ignore_codimensional_vertices) {
            hash_grid.addVerticesFromEdges(V, E, inflation_radius);
        } else {
            hash_grid.addVertices(V, inflation_radius);
        }
        hash_grid.addEdges(V, E, inflation_radius);
        if (dim == 3) {
            // These are not needed for 2D
            hash_grid.addFaces(V, F, inflation_radius);
        }

        if (dim == 2) {
            // This is not needed for 3D
            hash_grid.getVertexEdgePairs(
                E, candidates.ev_candidates, can_collide);
        } else {
            // These are not needed for 2D
            hash_grid.getEdgeEdgePairs(
                E, candidates.ee_candidates, can_collide);
            hash_grid.getFaceVertexPairs(
                F, candidates.fv_candidates, can_collide);
        }
    } break;
    case BroadPhaseMethod::SPATIAL_HASH: {
        // TODO: Use ignore_codimensional_vertices
        // TODO: Use can_collide
        SpatialHash sh(V, E, F, inflation_radius);
        sh.queryMeshForCandidates(
            V, E, F, candidates,
            /*queryEV=*/dim == 2,
            /*queryEE=*/dim == 3,
            /*queryFV=*/dim == 3);
    } break;
    }
}

void construct_collision_candidates(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    Candidates& candidates,
    double inflation_radius,
    const BroadPhaseMethod& method,
    bool ignore_codimensional_vertices,
    const std::function<bool(size_t, size_t)>& can_collide)
{
    int dim = V0.cols();
    assert(V1.cols() == dim);

    candidates.clear();

    switch (method) {
    case BroadPhaseMethod::BRUTE_FORCE:
        detect_collision_candidates_brute_force(
            V0, V1, E, F, candidates,
            /*detect_edge_vertex=*/dim == 2,
            /*detect_edge_edge=*/dim == 3,
            /*detect_face_vertex=*/dim == 3,
            /*perform_aabb_check=*/true,
            /*aabb_inflation_radius=*/inflation_radius,
            ignore_codimensional_vertices, can_collide);
        break;
    case BroadPhaseMethod::HASH_GRID: {
        HashGrid hash_grid;
        hash_grid.resize(V0, V1, E, inflation_radius);

        // Assumes the edges connect to all boundary vertices
        if (ignore_codimensional_vertices) {
            hash_grid.addVerticesFromEdges(V0, V1, E, inflation_radius);
        } else {
            hash_grid.addVertices(V0, V1, inflation_radius);
        }
        hash_grid.addEdges(V0, V1, E, inflation_radius);
        if (dim == 3) {
            // These are not needed for 2D
            hash_grid.addFaces(V0, V1, F, inflation_radius);
        }

        if (dim == 2) {
            // This is not needed for 3D
            hash_grid.getVertexEdgePairs(
                E, candidates.ev_candidates, can_collide);
        } else {
            // These are not needed for 2D
            hash_grid.getEdgeEdgePairs(
                E, candidates.ee_candidates, can_collide);
            hash_grid.getFaceVertexPairs(
                F, candidates.fv_candidates, can_collide);
        }
    } break;
    case BroadPhaseMethod::SPATIAL_HASH: {
        // TODO: Use ignore_codimensional_vertices
        // TODO: Use can_collide
        SpatialHash sh(V0, V1, E, F, inflation_radius);
        sh.queryMeshForCandidates(
            V0, V1, E, F, candidates,
            /*queryEV=*/dim == 2,
            /*queryEE=*/dim == 3,
            /*queryFV=*/dim == 3);
    } break;
    }
}

} // namespace ipc
