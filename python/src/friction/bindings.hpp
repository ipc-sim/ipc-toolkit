#pragma once

#include <friction/collisions/bindings.hpp>

#include <pybind11/pybind11.h>
namespace py = pybind11;

void define_closest_point(py::module_& m);
void define_friction_collisions(py::module_& m);
void define_normal_force_magnitude(py::module_& m);
void define_relative_velocity(py::module_& m);
void define_smooth_friction_mollifier(py::module_& m);
void define_tangent_basis(py::module_& m);