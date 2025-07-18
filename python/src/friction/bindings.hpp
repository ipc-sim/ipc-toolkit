#pragma once

#include <pybind11/pybind11.h>

void define_smooth_friction_mollifier(py::module_& m);
void define_smooth_mu(py::module_& m);