#pragma once

#include <pybind11/pybind11.h>

void define_distance_type(py::module_& m);
void define_edge_edge_mollifier(py::module_& m);
void define_edge_edge_distance(py::module_& m);
void define_line_line_distance(py::module_& m);
void define_point_edge_distance(py::module_& m);
void define_point_line_distance(py::module_& m);
void define_point_point_distance(py::module_& m);
void define_point_plane_distance(py::module_& m);
void define_point_triangle_distance(py::module_& m);
void define_signed_distance(py::module_& m);