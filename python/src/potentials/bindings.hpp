#pragma once

#include <pybind11/pybind11.h>
namespace py = pybind11;

void define_barrier_potential(py::module_& m);
void define_smooth_potential(py::module_& m);
void define_friction_potential(py::module_& m);
void define_normal_adhesion_potential(py::module_& m);
void define_normal_potential(py::module_& m);
void define_tangential_adhesion_potential(py::module_& m);
void define_tangential_potential(py::module_& m);
