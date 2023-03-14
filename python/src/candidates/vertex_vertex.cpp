#include <common.hpp>

#include <ipc/candidates/vertex_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_vertex_vertex_candidate(py::module_& m)
{
    py::class_<VertexVertexCandidate, CollisionStencil>(
        m, "VertexVertexCandidate")
        .def(
            py::init<long, long>(), "", py::arg("vertex0_id"),
            py::arg("vertex1_id"))
        .def(
            "__str__",
            [](const VertexVertexCandidate& ev) {
                return fmt::format(
                    "[{:d}, {:d}]", ev.vertex0_id, ev.vertex1_id);
            })
        .def(
            "__repr__",
            [](const VertexVertexCandidate& ev) {
                return fmt::format(
                    "VertexVertexCandidate({:d}, {:d})", ev.vertex0_id,
                    ev.vertex1_id);
            })
        .def("num_vertices", &VertexVertexCandidate::num_vertices, "")
        .def(
            "vertex_ids", &VertexVertexCandidate::vertex_ids,
            R"ipc_Qu8mg5v7(
            Get the indices of the vertices

            Parameters:
                edges: edge matrix of mesh
                faces: face matrix of mesh

            Returns:
                List of vertex indices
            )ipc_Qu8mg5v7",
            py::arg("edges"), py::arg("faces"))
        .def("__eq__", &VertexVertexCandidate::operator==, "", py::arg("other"))
        .def("__ne__", &VertexVertexCandidate::operator!=, "", py::arg("other"))
        .def(
            "__lt__", &VertexVertexCandidate::operator<,
            "Compare EdgeVertexCandidates for sorting.", py::arg("other"))
        .def_readwrite(
            "vertex0_id", &VertexVertexCandidate::vertex0_id,
            "ID of the first vertex")
        .def_readwrite(
            "vertex1_id", &VertexVertexCandidate::vertex1_id,
            "ID of the second vertex");
}
