#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

void define_adaptive_stiffness(py::module_& m);
void define_barrier(py::module_& m);