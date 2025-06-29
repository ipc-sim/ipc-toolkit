#include <common.hpp>

#include <ipc/candidates/edge_edge.hpp>

using namespace ipc;

void define_edge_edge_candidate(py::module_& m)
{
    py::class_<EdgeEdgeCandidate, CollisionStencil>(m, "EdgeEdgeCandidate")
        .def(py::init<index_t, index_t>(), "edge0_id"_a, "edge1_id"_a)
        .def(
            py::init([](std::tuple<index_t, index_t> edge_ids) {
                return std::make_unique<EdgeEdgeCandidate>(
                    std::get<0>(edge_ids), std::get<1>(edge_ids));
            }),
            "edge_ids"_a)
        .def("known_dtype", &EdgeEdgeCandidate::known_dtype)
        .def(
            "__str__",
            [](const EdgeEdgeCandidate& ee) {
                return fmt::format("[{:d}, {:d}]", ee.edge0_id, ee.edge1_id);
            })
        .def(
            "__repr__",
            [](const EdgeEdgeCandidate& ee) {
                return fmt::format(
                    "EdgeEdgeCandidate({:d}, {:d})", ee.edge0_id, ee.edge1_id);
            })
        .def("__eq__", &EdgeEdgeCandidate::operator==, "other"_a)
        .def("__ne__", &EdgeEdgeCandidate::operator!=, "other"_a)
        .def(
            "__lt__", &EdgeEdgeCandidate::operator<,
            "Compare EdgeEdgeCandidates for sorting.", "other"_a)
        .def_readwrite(
            "edge0_id", &EdgeEdgeCandidate::edge0_id, "ID of the first edge.")
        .def_readwrite(
            "edge1_id", &EdgeEdgeCandidate::edge1_id, "ID of the second edge.");

    py::implicitly_convertible<
        std::tuple<index_t, index_t>, EdgeEdgeCandidate>();
}
