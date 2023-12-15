#pragma once

#include <barrier/bindings.hpp>
#include <broad_phase/bindings.hpp>
#include <candidates/bindings.hpp>
#include <ccd/bindings.hpp>
#include <collisions/bindings.hpp>
#include <distance/bindings.hpp>
#include <friction/bindings.hpp>
#include <implicits/bindings.hpp>
#include <potentials/bindings.hpp>
#include <utils/bindings.hpp>

#include <pybind11/pybind11.h>

namespace py = pybind11;

void define_collision_mesh(py::module_& m);
void define_ipc(py::module_& m);