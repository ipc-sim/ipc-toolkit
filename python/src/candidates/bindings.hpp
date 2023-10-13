#pragma once

#include <pybind11/pybind11.h>

// candidates
void define_candidates(py::module_& m);
void define_collision_stencil(py::module_& m);
void define_continuous_collision_candidate(py::module_& m);
void define_edge_edge_candidate(py::module_& m);
void define_edge_face_candidate(py::module_& m);
void define_edge_vertex_candidate(py::module_& m);
void define_face_vertex_candidate(py::module_& m);
void define_vertex_vertex_candidate(py::module_& m);