#include <common.hpp>

#include <ipc/candidates/vertex_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_vertex_vertex_candidate(py::module_& m)
{
    py::class_<VertexVertexCandidate, CollisionStencil>(
        m, "VertexVertexCandidate")
        .def(
            py::init<index_t, index_t>(), py::arg("vertex0_id"),
            py::arg("vertex1_id"))
        .def(
            py::init([](std::tuple<index_t, index_t> vertex_ids) {
                return std::make_unique<VertexVertexCandidate>(
                    std::get<0>(vertex_ids), std::get<1>(vertex_ids));
            }),
            py::arg("vertex_ids"))
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
        .def("__eq__", &VertexVertexCandidate::operator==, py::arg("other"))
        .def("__ne__", &VertexVertexCandidate::operator!=, py::arg("other"))
        .def(
            "__lt__", &VertexVertexCandidate::operator<,
            "Compare EdgeVertexCandidates for sorting.", py::arg("other"))
        .def_readwrite(
            "vertex0_id", &VertexVertexCandidate::vertex0_id,
            "ID of the first vertex")
        .def_readwrite(
            "vertex1_id", &VertexVertexCandidate::vertex1_id,
            "ID of the second vertex");

    py::implicitly_convertible<
        std::tuple<index_t, index_t>, VertexVertexCandidate>();
}
