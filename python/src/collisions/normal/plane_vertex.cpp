#include <common.hpp>

#include <ipc/collisions/normal/plane_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_plane_vertex_normal_collision(py::module_& m)
{
    py::class_<PlaneVertexNormalCollision, NormalCollision>(
        m, "PlaneVertexNormalCollision")
        .def(
            py::init<const VectorMax3d&, const VectorMax3d&, const long>(),
            py::arg("plane_origin"), py::arg("plane_normal"),
            py::arg("vertex_id"))
        .def_readwrite(
            "plane_origin", &PlaneVertexNormalCollision::plane_origin,
            "The plane's origin.")
        .def_readwrite(
            "plane_normal", &PlaneVertexNormalCollision::plane_normal,
            "The plane's normal.")
        .def_readwrite(
            "vertex_id", &PlaneVertexNormalCollision::vertex_id,
            "The vertex's id.");
}
