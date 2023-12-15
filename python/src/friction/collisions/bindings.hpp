#pragma once

#include <pybind11/pybind11.h>
namespace py = pybind11;

void define_edge_edge_friction_collision(py::module_& m);
void define_edge_vertex_friction_collision(py::module_& m);
void define_face_vertex_friction_collision(py::module_& m);
void define_friction_collision(py::module_& m);
void define_vertex_vertex_friction_collision(py::module_& m);