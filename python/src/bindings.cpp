#include <common.hpp>
#include <bindings.hpp>

PYBIND11_MODULE(ipctk, m)
{
    // py::options options;
    // options.disable_function_signatures();

    m.doc() = "IPC Toolkit";

    // define these early because they are used in other definitions
    define_eigen_ext(m);
    define_narrow_phase_ccd(m);
    define_tight_inclusion_ccd(m);

    // adhesion
    define_adhesion(m);

    // barrier
    define_barrier(m);
    define_barrier_force_magnitude(m);
    define_adaptive_stiffness(m);

    // broad_phase
    define_aabb(m);
    define_broad_phase(m);
    define_brute_force(m);
    define_bvh(m);
    define_hash_grid(m);
    define_lbvh(m);
    define_spatial_hash(m);
    define_sweep_and_prune(m);
    define_sweep_and_tiniest_queue(m);
    define_voxel_size_heuristic(m);

    // candidates
    define_candidates(m);
    define_collision_stencil(m);
    define_edge_edge_candidate(m);
    define_edge_face_candidate(m);
    define_edge_vertex_candidate(m);
    define_face_face_candidate(m);
    define_face_vertex_candidate(m);
    define_vertex_vertex_candidate(m);

    // ccd
    define_ccd_aabb(m);
    define_check_initial_distance(m);
    define_inexact_point_edge(m);
    define_point_static_plane(m);
    define_inexact_ccd(m);
    define_additive_ccd(m);
    define_nonlinear_ccd(m);

    // collisions/normal
    define_distance_type(m); // define early because it is used next
    define_normal_collision(m);
    define_normal_collisions(m);
    define_edge_edge_normal_collision(m);
    define_edge_vertex_normal_collision(m);
    define_face_vertex_normal_collision(m);
    define_plane_vertex_normal_collision(m);
    define_vertex_vertex_normal_collision(m);

    // tangent
    define_closest_point(m);
    define_relative_velocity(m);
    define_tangent_basis(m);

    // collisions/tangential
    define_tangential_collision(m);
    define_tangential_collisions(m);
    define_edge_edge_tangential_collision(m);
    define_edge_vertex_tangential_collision(m);
    define_face_vertex_tangential_collision(m);
    define_vertex_vertex_tangential_collision(m);

    // distance
    define_edge_edge_mollifier(m);
    define_edge_edge_distance(m);
    define_line_line_distance(m);
    define_point_edge_distance(m);
    define_point_line_distance(m);
    define_point_point_distance(m);
    define_point_plane_distance(m);
    define_point_triangle_distance(m);
    define_signed_distance(m);

    // friction
    define_smooth_friction_mollifier(m);
    define_smooth_mu(m);

    define_smooth_potential(m);

    // geometry
    define_angle(m);
    define_area(m);
    define_intersection(m);
    define_normal(m);

    // implicits
    define_plane_implicit(m);

    // math
    define_interval(m);

    // ogc
    py::module_ ogc =
        m.def_submodule("ogc", "Offset Geometric Contact (OGC) helpers");
    define_feasible_region(ogc);
    define_trust_region(ogc);

    // potentials
    define_normal_potential(m); // define early because it is used next
    define_barrier_potential(m);
    define_normal_adhesion_potential(m);
    define_tangential_potential(m);
    define_friction_potential(m);
    define_tangential_adhesion_potential(m);

    // utils
    define_logger(m);
    define_thread_limiter(m);
    define_vertex_to_min_edge(m);
    define_world_bbox_diagonal_length(m);

    // root
    define_collision_mesh(m);
    define_ipc(m);
}
