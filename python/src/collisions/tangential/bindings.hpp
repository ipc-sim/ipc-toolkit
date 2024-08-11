#pragma once

#include <pybind11/pybind11.h>
namespace py = pybind11;

void define_tangential_collision(py::module_& m);
void define_tangential_collisions(py::module_& m);
void define_edge_edge_tangential_collision(py::module_& m);
void define_edge_vertex_tangential_collision(py::module_& m);
void define_face_vertex_tangential_collision(py::module_& m);
void define_vertex_vertex_tangential_collision(py::module_& m);