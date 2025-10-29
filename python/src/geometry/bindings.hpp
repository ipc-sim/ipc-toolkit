#pragma once

#include <pybind11/pybind11.h>

void define_angle(py::module_& m);
void define_area(py::module_& m);
void define_normal(py::module_& m);
void define_intersection(py::module_& m);