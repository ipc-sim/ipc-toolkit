#include <common.hpp>

#include <ipc/candidates/edge_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_edge_vertex_candidate(py::module_& m)
{
    py::class_<EdgeVertexCandidate, CollisionStencil>(m, "EdgeVertexCandidate")
        .def(py::init<long, long>(), py::arg("edge_id"), py::arg("vertex_id"))
        .def(
            py::init([](std::tuple<long, long> edge_and_vertex_id) {
                return std::make_unique<EdgeVertexCandidate>(
                    std::get<0>(edge_and_vertex_id),
                    std::get<1>(edge_and_vertex_id));
            }),
            py::arg("edge_and_vertex_id"))
        .def("known_dtype", &EdgeVertexCandidate::known_dtype)
        .def(
            "__str__",
            [](const EdgeVertexCandidate& ev) {
                return fmt::format("[{:d}, {:d}]", ev.edge_id, ev.vertex_id);
            })
        .def(
            "__repr__",
            [](const EdgeVertexCandidate& ev) {
                return fmt::format(
                    "EdgeVertexCandidate({:d}, {:d})", ev.edge_id,
                    ev.vertex_id);
            })
        .def("__eq__", &EdgeVertexCandidate::operator==, py::arg("other"))
        .def("__ne__", &EdgeVertexCandidate::operator!=, py::arg("other"))
        .def(
            "__lt__", &EdgeVertexCandidate::operator<,
            "Compare EdgeVertexCandidates for sorting.", py::arg("other"))
        .def_readwrite(
            "edge_id", &EdgeVertexCandidate::edge_id, "ID of the edge")
        .def_readwrite(
            "vertex_id", &EdgeVertexCandidate::vertex_id, "ID of the vertex");

    py::implicitly_convertible<std::tuple<long, long>, EdgeVertexCandidate>();
}
