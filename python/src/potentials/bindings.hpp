#pragma once

#include <pybind11/pybind11.h>
namespace py = pybind11;

void define_barrier_potential(py::module_& m);
void define_friction_potential(py::module_& m);