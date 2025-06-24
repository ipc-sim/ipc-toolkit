#include <common.hpp>

#include <ipc/candidates/edge_vertex.hpp>

using namespace ipc;

void define_edge_vertex_candidate(py::module_& m)
{
    py::class_<EdgeVertexCandidate, CollisionStencil>(m, "EdgeVertexCandidate")
        .def(py::init<index_t, index_t>(), "edge_id"_a, "vertex_id"_a)
        .def(
            py::init([](std::tuple<index_t, index_t> edge_and_vertex_id) {
                return std::make_unique<EdgeVertexCandidate>(
                    std::get<0>(edge_and_vertex_id),
                    std::get<1>(edge_and_vertex_id));
            }),
            "edge_and_vertex_id"_a)
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
        .def("__eq__", &EdgeVertexCandidate::operator==, "other"_a)
        .def("__ne__", &EdgeVertexCandidate::operator!=, "other"_a)
        .def(
            "__lt__", &EdgeVertexCandidate::operator<,
            "Compare EdgeVertexCandidates for sorting.", "other"_a)
        .def_readwrite(
            "edge_id", &EdgeVertexCandidate::edge_id, "ID of the edge")
        .def_readwrite(
            "vertex_id", &EdgeVertexCandidate::vertex_id, "ID of the vertex");

    py::implicitly_convertible<
        std::tuple<index_t, index_t>, EdgeVertexCandidate>();
}
