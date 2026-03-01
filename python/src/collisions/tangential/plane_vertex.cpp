#include <ipc/collisions/tangential/plane_vertex.hpp>

#include <common.hpp>

using namespace ipc;

void define_plane_vertex_tangential_collision(py::module_& m)
{
    py::class_<
        PlaneVertexTangentialCollision, PlaneVertexCandidate,
        TangentialCollision>(m, "PlaneVertexTangentialCollision")
        .def(py::init<const PlaneVertexNormalCollision&>(), "collision"_a)
        .def(
            py::init<
                const PlaneVertexNormalCollision&,
                Eigen::ConstRef<VectorMax12d>, const NormalPotential&>(),
            "collision"_a, "positions"_a, "normal_potential"_a);
}