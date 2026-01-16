#pragma once

#include <pybind11/pybind11.h>

void define_feasible_region(py::module_& m);
void define_trust_region(py::module_& m);