#include <common.hpp>

#include <ipc/collisions/normal/vertex_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_vertex_vertex_normal_collision(py::module_& m)
{
    py::class_<
        VertexVertexNormalCollision, VertexVertexCandidate, NormalCollision>(
        m, "VertexVertexNormalCollision")
        .def(py::init<const VertexVertexCandidate&>(), py::arg("vv_candidate"))
        .def(
            py::init<long, long>(), "", py::arg("vertex0_id"),
            py::arg("vertex1_id"));
}
