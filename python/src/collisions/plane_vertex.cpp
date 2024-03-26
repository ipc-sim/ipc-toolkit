#include <common.hpp>

#include <ipc/collisions/plane_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_plane_vertex_collision(py::module_& m)
{
    py::class_<PlaneVertexCollision, Collision>(m, "PlaneVertexCollision")
        .def(
            py::init<const VectorMax3d&, const VectorMax3d&, const long>(),
            py::arg("plane_origin"), py::arg("plane_normal"),
            py::arg("vertex_id"))
        .def_readwrite("plane_origin", &PlaneVertexCollision::plane_origin)
        .def_readwrite("plane_normal", &PlaneVertexCollision::plane_normal)
        .def_readwrite("vertex_id", &PlaneVertexCollision::vertex_id);
}
