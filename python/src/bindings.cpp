#include <pybind11/pybind11.h>

namespace py = pybind11;

void define_barrier_functions(py::module_& m);

void define_broad_phase_functions(py::module_& m);

void define_collision_candidate_classes(py::module_& m);
void define_collision_constraint_classes(py::module_& m);

void define_distance_type_functions(py::module_& m);
void define_edge_edge_mollifier_functions(py::module_& m);
void define_edge_edge_distance_functions(py::module_& m);
void define_line_line_distance_functions(py::module_& m);
void define_point_edge_distance_functions(py::module_& m);
void define_point_line_distance_functions(py::module_& m);
void define_point_point_distance_functions(py::module_& m);
void define_point_plane_distance_functions(py::module_& m);
void define_point_triangle_distance_functions(py::module_& m);

void define_ipc_functions(py::module_& m);

void define_collision_mesh_class(py::module_& m);

void define_logger_functions(py::module_& m);
void define_thread_limiter_functions(py::module_& m);

PYBIND11_MODULE(ipctk, m)
{
    m.doc() = "IPC Toolkit";

    define_barrier_functions(m);

    define_broad_phase_functions(m);

    define_collision_candidate_classes(m);
    define_collision_constraint_classes(m);

    define_distance_type_functions(m);
    define_edge_edge_mollifier_functions(m);
    define_edge_edge_distance_functions(m);
    define_line_line_distance_functions(m);
    define_point_edge_distance_functions(m);
    define_point_line_distance_functions(m);
    define_point_point_distance_functions(m);
    define_point_plane_distance_functions(m);
    define_point_triangle_distance_functions(m);

    define_ipc_functions(m);

    define_collision_mesh_class(m);

    define_logger_functions(m);
    define_thread_limiter_functions(m);
}
