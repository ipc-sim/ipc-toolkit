#pragma once

#include <pybind11/pybind11.h>
namespace py = pybind11;

void define_ccd_aabb(py::module_& m);
void define_ccd(py::module_& m);
void define_additive_ccd(py::module_& m);
void define_inexact_point_edge(py::module_& m);
void define_nonlinear_ccd(py::module_& m);
void define_point_static_plane(py::module_& m);
void define_tight_inclusion_ccd(py::module_& m);