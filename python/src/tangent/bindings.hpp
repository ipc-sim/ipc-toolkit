#pragma once

#include <pybind11/pybind11.h>
namespace py = pybind11;

void define_closest_point(py::module_& m);
void define_relative_velocity(py::module_& m);
void define_tangent_basis(py::module_& m);