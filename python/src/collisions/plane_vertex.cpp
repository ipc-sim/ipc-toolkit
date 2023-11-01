#include <common.hpp>

#include <ipc/collisions/plane_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_plane_vertex_constraint(py::module_& m)
{
    py::class_<PlaneVertexConstraint, CollisionConstraint>(
        m, "PlaneVertexConstraint")
        .def(
            py::init<const VectorMax3d&, const VectorMax3d&, const long>(),
            py::arg("plane_origin"), py::arg("plane_normal"),
            py::arg("vertex_id"))
        .def_readwrite("plane_origin", &PlaneVertexConstraint::plane_origin)
        .def_readwrite("plane_normal", &PlaneVertexConstraint::plane_normal)
        .def_readwrite("vertex_id", &PlaneVertexConstraint::vertex_id);
}
