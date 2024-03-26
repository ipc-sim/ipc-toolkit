#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

void define_aabb(py::module_& m);
void define_broad_phase(py::module_& m);
void define_brute_force(py::module_& m);
void define_bvh(py::module_& m);
void define_hash_grid(py::module_& m);
void define_spatial_hash(py::module_& m);
void define_sweep_and_prune(py::module_& m);
void define_sweep_and_tiniest_queue(py::module_& m);
void define_voxel_size_heuristic(py::module_& m);
