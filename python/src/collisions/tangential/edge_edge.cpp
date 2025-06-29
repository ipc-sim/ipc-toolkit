#include <common.hpp>

#include <ipc/collisions/tangential/edge_edge.hpp>

using namespace ipc;

void define_edge_edge_tangential_collision(py::module_& m)
{
    py::class_<
        EdgeEdgeTangentialCollision, EdgeEdgeCandidate, TangentialCollision>(
        m, "EdgeEdgeTangentialCollision")
        .def(py::init<const EdgeEdgeNormalCollision&>(), "collision"_a)
        .def(
            py::init<
                const EdgeEdgeNormalCollision&, Eigen::ConstRef<VectorMax12d>,
                const NormalPotential&, const double>(),
            "collision"_a, "positions"_a, "normal_potential"_a,
            "normal_stiffness"_a);
}
