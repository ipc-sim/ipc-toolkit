#include "../common.hpp"

#include <ipc/candidates/vertex_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_vertex_vertex_candidate(py::module_& m)
{
    py::class_<VertexVertexCandidate>(m, "VertexVertexCandidate")
        .def(
            py::init<long, long>(), "", py::arg("vertex0_index"),
            py::arg("vertex1_index"))
        .def(
            "__str__",
            [](const VertexVertexCandidate& ev) {
                return fmt::format(
                    "[{:d}, {:d}]", ev.vertex0_index, ev.vertex1_index);
            })
        .def(
            "__repr__",
            [](const VertexVertexCandidate& ev) {
                return fmt::format(
                    "VertexVertexCandidate({:d}, {:d})", ev.vertex0_index,
                    ev.vertex1_index);
            })
        .def("num_vertices", &VertexVertexCandidate::num_vertices, "")
        .def(
            "vertex_indices", &VertexVertexCandidate::vertex_indices,
            R"ipc_Qu8mg5v7(
            Get the indices of the vertices

            Parameters:
                E: edge matrix of mesh
                F: face matrix of mesh

            Returns:
                List of vertex indices
            )ipc_Qu8mg5v7",
            py::arg("E"), py::arg("F"))
        .def(
            "compute_distance", &VertexVertexCandidate::compute_distance, "",
            py::arg("V"), py::arg("E"), py::arg("F"))
        .def(
            "compute_distance_gradient",
            &VertexVertexCandidate::compute_distance_gradient, "", py::arg("V"),
            py::arg("E"), py::arg("F"))
        .def(
            "compute_distance_hessian",
            &VertexVertexCandidate::compute_distance_hessian, "", py::arg("V"),
            py::arg("E"), py::arg("F"))
        .def("__eq__", &VertexVertexCandidate::operator==, "", py::arg("other"))
        .def("__ne__", &VertexVertexCandidate::operator!=, "", py::arg("other"))
        .def(
            "__lt__", &VertexVertexCandidate::operator<,
            "Compare EdgeVertexCandidates for sorting.", py::arg("other"))
        .def_readwrite(
            "vertex0_index", &VertexVertexCandidate::vertex0_index, "")
        .def_readwrite(
            "vertex1_index", &VertexVertexCandidate::vertex1_index, "");
}
