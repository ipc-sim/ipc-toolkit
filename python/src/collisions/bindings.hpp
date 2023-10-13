#pragma once

#include <pybind11/pybind11.h>
namespace py = pybind11;

void define_collision_constraint(py::module_& m);
void define_collision_constraints(py::module_& m);
void define_edge_edge_constraint(py::module_& m);
void define_edge_vertex_constraint(py::module_& m);
void define_face_vertex_constraint(py::module_& m);
void define_plane_vertex_constraint(py::module_& m);
void define_vertex_vertex_constraint(py::module_& m);