#pragma once

#include <pybind11/pybind11.h>
namespace py = pybind11;

void define_area_gradient(py::module_& m);
void define_eigen_ext(py::module_& m);
void define_interval(py::module_& m);
void define_intersection(py::module_& m);
void define_logger(py::module_& m);
void define_thread_limiter(py::module_& m);
void define_vertex_to_min_edge(py::module_& m);
void define_world_bbox_diagonal_length(py::module_& m);