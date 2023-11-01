#include <common.hpp>
#include <bindings.hpp>

namespace py = pybind11;

PYBIND11_MODULE(ipctk, m)
{
    // py::options options;
    // options.disable_function_signatures();

    m.doc() = "IPC Toolkit";

    // barrier
    define_adaptive_stiffness(m);
    define_barrier(m);

    // broad_phase
    define_aabb(m);
    define_broad_phase(m);
    define_brute_force(m);
    define_bvh(m);
    define_hash_grid(m);
    define_spatial_hash(m);
    define_sweep_and_tiniest_queue(m);
    define_voxel_size_heuristic(m);

    // candidates
    define_candidates(m);
    define_collision_stencil(m);
    define_continuous_collision_candidate(m);
    define_edge_edge_candidate(m);
    define_edge_face_candidate(m);
    define_edge_vertex_candidate(m);
    define_face_vertex_candidate(m);
    define_vertex_vertex_candidate(m);

    // ccd
    define_ccd_aabb(m);
    define_ccd(m);
    define_additive_ccd(m);
    define_inexact_point_edge(m);
    define_nonlinear_ccd(m);
    define_point_static_plane(m);

    // collisions
    define_distance_type(m); // define early because it is used next
    define_collision_constraint(m);
    define_collision_constraints(m);
    define_edge_edge_constraint(m);
    define_edge_vertex_constraint(m);
    define_face_vertex_constraint(m);
    define_plane_vertex_constraint(m);
    define_vertex_vertex_constraint(m);

    // distance
    define_edge_edge_mollifier(m);
    define_edge_edge_distance(m);
    define_line_line_distance(m);
    define_point_edge_distance(m);
    define_point_line_distance(m);
    define_point_point_distance(m);
    define_point_plane_distance(m);
    define_point_triangle_distance(m);

    // friction
    define_closest_point(m);
    define_friction_constraints(m);
    define_normal_force_magnitude(m);
    define_relative_velocity(m);
    define_smooth_friction_mollifier(m);
    define_tangent_basis(m);

    // friction/constraints
    // NOTE: this has to be defined before the other friction constraints
    define_friction_constraint(m);
    define_edge_edge_friction_constraint(m);
    define_edge_vertex_friction_constraint(m);
    define_face_vertex_friction_constraint(m);
    define_vertex_vertex_friction_constraint(m);

    // implicits
    define_plane_implicit(m);

    // utils
    define_area_gradient(m);
    define_eigen_ext(m);
    define_interval(m);
    define_intersection(m);
    define_logger(m);
    define_thread_limiter(m);
    define_vertex_to_min_edge(m);
    define_world_bbox_diagonal_length(m);

    // root
    define_collision_mesh(m);
    define_ipc(m);
}
