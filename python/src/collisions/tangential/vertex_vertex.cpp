#include <common.hpp>

#include <ipc/collisions/tangential/vertex_vertex.hpp>

using namespace ipc;

void define_vertex_vertex_tangential_collision(py::module_& m)
{
    py::class_<
        VertexVertexTangentialCollision, VertexVertexCandidate,
        TangentialCollision>(m, "VertexVertexTangentialCollision")
        .def(py::init<const VertexVertexNormalCollision&>(), "collision"_a)
        .def(
            py::init<
                const VertexVertexNormalCollision&,
                Eigen::ConstRef<VectorMax12d>, const NormalPotential&,
                const double>(),
            "collision"_a, "positions"_a, "normal_potential"_a,
            "normal_stiffness"_a);
}
