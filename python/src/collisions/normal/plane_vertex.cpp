#include <Eigen/Geometry>
#include <common.hpp>

#include <ipc/collisions/normal/plane_vertex.hpp>

using namespace ipc;

void define_plane_vertex_normal_collision(py::module_& m)
{
    py::class_<
        PlaneVertexNormalCollision, PlaneVertexCandidate, NormalCollision>(
        m, "PlaneVertexNormalCollision")
        .def(
            py::init<Eigen::Hyperplane<double, 3>, index_t>(), "plane"_a,
            "vertex_id"_a);
}
