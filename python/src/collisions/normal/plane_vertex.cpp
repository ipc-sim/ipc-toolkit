#include <common.hpp>

#include <ipc/collisions/normal/plane_vertex.hpp>

using namespace ipc;

void define_plane_vertex_normal_collision(py::module_& m)
{
    py::class_<PlaneVertexNormalCollision, NormalCollision>(
        m, "PlaneVertexNormalCollision")
        .def(
            py::init<
                Eigen::ConstRef<VectorMax3d>, Eigen::ConstRef<VectorMax3d>,
                const index_t>(),
            "plane_origin"_a, "plane_normal"_a, "vertex_id"_a)
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
