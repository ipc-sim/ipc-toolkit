#include <common.hpp>

#include <ipc/collisions/vertex_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_vertex_vertex_constraint(py::module_& m)
{
    py::class_<
        VertexVertexConstraint, VertexVertexCandidate, CollisionConstraint>(
        m, "VertexVertexConstraint")
        .def(
            py::init<long, long>(), "", py::arg("vertex0_id"),
            py::arg("vertex1_id"))
        .def(py::init<VertexVertexCandidate>(), "", py::arg("vv_candidate"));
}
