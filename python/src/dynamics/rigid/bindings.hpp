#pragma once

#include <pybind11/pybind11.h>
namespace py = pybind11;

void define_rigid_simulator(py::module_& m);