#include <common.hpp>

#include <ipc/collisions/tangential/edge_vertex.hpp>

using namespace ipc;

void define_edge_vertex_tangential_collision(py::module_& m)
{
    py::class_<
        EdgeVertexTangentialCollision, EdgeVertexCandidate,
        TangentialCollision>(m, "EdgeVertexTangentialCollision")
        .def(py::init<const EdgeVertexNormalCollision&>(), "collision"_a)
        .def(
            py::init<
                const EdgeVertexNormalCollision&, Eigen::ConstRef<VectorMax12d>,
                const NormalPotential&, const double>(),
            "collision"_a, "positions"_a, "normal_potential"_a,
            "normal_stiffness"_a);
}
