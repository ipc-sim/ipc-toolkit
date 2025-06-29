#include <common.hpp>

#include <ipc/candidates/vertex_vertex.hpp>

using namespace ipc;

void define_vertex_vertex_candidate(py::module_& m)
{
    py::class_<VertexVertexCandidate, CollisionStencil>(
        m, "VertexVertexCandidate")
        .def(py::init<index_t, index_t>(), "vertex0_id"_a, "vertex1_id"_a)
        .def(
            py::init([](std::tuple<index_t, index_t> vertex_ids) {
                return std::make_unique<VertexVertexCandidate>(
                    std::get<0>(vertex_ids), std::get<1>(vertex_ids));
            }),
            "vertex_ids"_a)
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
        .def("__eq__", &VertexVertexCandidate::operator==, "other"_a)
        .def("__ne__", &VertexVertexCandidate::operator!=, "other"_a)
        .def(
            "__lt__", &VertexVertexCandidate::operator<,
            "Compare EdgeVertexCandidates for sorting.", "other"_a)
        .def_readwrite(
            "vertex0_id", &VertexVertexCandidate::vertex0_id,
            "ID of the first vertex")
        .def_readwrite(
            "vertex1_id", &VertexVertexCandidate::vertex1_id,
            "ID of the second vertex");

    py::implicitly_convertible<
        std::tuple<index_t, index_t>, VertexVertexCandidate>();
}
