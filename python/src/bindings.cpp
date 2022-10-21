#include "common.hpp"

namespace py = pybind11;

// barrier
void define_adaptive_stiffness(py::module_& m);
void define_barrier(py::module_& m);

// broad_phase
void define_aabb(py::module_& m);
void define_broad_phase(py::module_& m);
void define_brute_force(py::module_& m);
void define_hash_grid(py::module_& m);
void define_spatial_hash(py::module_& m);
void define_sweep(py::module_& m);
void define_voxel_size_heuristic(py::module_& m);

// candidates
void define_candidates(py::module_& m);
void define_continuous_collision_candidate(py::module_& m);
void define_edge_edge_candidate(py::module_& m);
void define_edge_face_candidate(py::module_& m);
void define_edge_vertex_candidate(py::module_& m);
void define_face_vertex_candidate(py::module_& m);
void define_vertex_vertex_candidate(py::module_& m);

// ccd
void define_ccd_aabb(py::module_& m);
void define_ccd(py::module_& m);
void define_inexact_point_edge(py::module_& m);
void define_point_static_plane(py::module_& m);

// distance
void define_distance_type(py::module_& m);
void define_edge_edge_mollifier(py::module_& m);
void define_edge_edge_distance(py::module_& m);
void define_line_line_distance(py::module_& m);
void define_point_edge_distance(py::module_& m);
void define_point_line_distance(py::module_& m);
void define_point_point_distance(py::module_& m);
void define_point_plane_distance(py::module_& m);
void define_point_triangle_distance(py::module_& m);

// friction
void define_closest_point(py::module_& m);
void define_friction_constraint(py::module_& m);
void define_friction(py::module_& m);
void define_normal_force_magnitude(py::module_& m);
void define_relative_displacement(py::module_& m);
void define_smooth_friction_mollifier(py::module_& m);
void define_tangent_basis(py::module_& m);

// utils
void define_eigen_ext(py::module_& m);
void define_intersection(py::module_& m);
void define_logger(py::module_& m);
void define_thread_limiter(py::module_& m);
void define_world_bbox_diagonal_length(py::module_& m);

// root
void define_collision_constraint(py::module_& m);
void define_collision_mesh(py::module_& m);
void define_ipc(py::module_& m);

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
    define_hash_grid(m);
    define_spatial_hash(m);
    define_sweep(m);
    define_voxel_size_heuristic(m);

    // distance
    define_distance_type(m);
    define_edge_edge_mollifier(m);
    define_edge_edge_distance(m);
    define_line_line_distance(m);
    define_point_edge_distance(m);
    define_point_line_distance(m);
    define_point_point_distance(m);
    define_point_plane_distance(m);
    define_point_triangle_distance(m);

    // candidates
    define_candidates(m);
    define_continuous_collision_candidate(m);
    define_edge_edge_candidate(m);
    define_edge_face_candidate(m);
    define_edge_vertex_candidate(m);
    define_face_vertex_candidate(m);
    define_vertex_vertex_candidate(m);

    // ccd
    define_ccd_aabb(m);
    define_ccd(m);
    define_inexact_point_edge(m);
    define_point_static_plane(m);

    // friction
    define_closest_point(m);
    define_friction_constraint(m);
    define_friction(m);
    define_normal_force_magnitude(m);
    define_relative_displacement(m);
    define_smooth_friction_mollifier(m);
    define_tangent_basis(m);

    // utils
    define_eigen_ext(m);
    define_intersection(m);
    define_logger(m);
    define_thread_limiter(m);
    define_world_bbox_diagonal_length(m);

    // root
    define_collision_constraint(m);
    define_collision_mesh(m);
    define_ipc(m);
}
