#pragma once

#include <pybind11/pybind11.h>
namespace py = pybind11;

void define_normal_collision(py::module_& m);
void define_normal_collisions(py::module_& m);
void define_edge_edge_normal_collision(py::module_& m);
void define_edge_vertex_normal_collision(py::module_& m);
void define_face_vertex_normal_collision(py::module_& m);
void define_plane_vertex_normal_collision(py::module_& m);
void define_vertex_vertex_normal_collision(py::module_& m);