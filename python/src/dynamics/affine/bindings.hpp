#pragma once

#include <pybind11/pybind11.h>
namespace py = pybind11;

void define_affine_joints(py::module_& m);
