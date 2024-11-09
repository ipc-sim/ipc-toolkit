#pragma once

#include <pybind11/pybind11.h>
namespace py = pybind11;

void define_smooth_friction_mollifier(py::module_& m);