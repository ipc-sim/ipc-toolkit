#pragma once

#include <pybind11/pybind11.h>

void define_adaptive_stiffness(py::module_& m);
void define_barrier_force_magnitude(py::module_& m);
void define_barrier(py::module_& m);