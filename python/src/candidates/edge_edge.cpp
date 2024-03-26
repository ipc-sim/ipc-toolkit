#include <common.hpp>

#include <ipc/candidates/edge_edge.hpp>

namespace py = pybind11;
using namespace ipc;

void define_edge_edge_candidate(py::module_& m)
{
    py::class_<
        EdgeEdgeCandidate, CollisionStencil, ContinuousCollisionCandidate>(
        m, "EdgeEdgeCandidate")
        .def(py::init<long, long>(), py::arg("edge0_id"), py::arg("edge1_id"))
        .def("known_dtype", &EdgeEdgeCandidate::known_dtype)
        .def(
            "__str__",
            [](const EdgeEdgeCandidate& ev) {
                return fmt::format("[{:d}, {:d}]", ev.edge0_id, ev.edge1_id);
            })
        .def(
            "__repr__",
            [](const EdgeEdgeCandidate& ev) {
                return fmt::format(
                    "EdgeEdgeCandidate({:d}, {:d})", ev.edge0_id, ev.edge1_id);
            })
        .def("__eq__", &EdgeEdgeCandidate::operator==, py::arg("other"))
        .def("__ne__", &EdgeEdgeCandidate::operator!=, py::arg("other"))
        .def(
            "__lt__", &EdgeEdgeCandidate::operator<,
            "Compare EdgeEdgeCandidates for sorting.", py::arg("other"))
        .def_readwrite(
            "edge0_id", &EdgeEdgeCandidate::edge0_id, "ID of the first edge.")
        .def_readwrite(
            "edge1_id", &EdgeEdgeCandidate::edge1_id, "ID of the second edge.");
}
