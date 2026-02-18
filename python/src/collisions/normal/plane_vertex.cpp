#include <common.hpp>

#include <ipc/collisions/normal/plane_vertex.hpp>

using namespace ipc;

void define_plane_vertex_normal_collision(py::module_& m)
{
    py::class_<
        PlaneVertexNormalCollision, PlaneVertexCandidate, NormalCollision>(
        m, "PlaneVertexNormalCollision")
        .def(
            py::init<
                Eigen::ConstRef<VectorMax3d>, Eigen::ConstRef<VectorMax3d>,
                index_t>(),
            "plane_origin"_a, "plane_normal"_a, "vertex_id"_a);
}
