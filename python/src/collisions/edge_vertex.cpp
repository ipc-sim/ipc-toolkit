#include <common.hpp>

#include <ipc/collisions/edge_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_edge_vertex_constraint(py::module_& m)
{
    py::class_<EdgeVertexConstraint, EdgeVertexCandidate, CollisionConstraint>(
        m, "EdgeVertexConstraint")
        .def(
            py::init<long, long>(), "", py::arg("edge_id"),
            py::arg("vertex_id"))
        .def(py::init<EdgeVertexCandidate>(), "", py::arg("ev_candidate"));
}
