#include <common.hpp>

#include <ipc/collisions/tangential/edge_edge.hpp>

namespace py = pybind11;
using namespace ipc;

void define_edge_edge_tangential_collision(py::module_& m)
{
    py::class_<
        EdgeEdgeTangentialCollision, EdgeEdgeCandidate, TangentialCollision>(
        m, "EdgeEdgeTangentialCollision")
        .def(py::init<const EdgeEdgeNormalCollision&>(), py::arg("collision"))
        .def(
            py::init<
                const EdgeEdgeNormalCollision&, const VectorMax12d&,
                const BarrierPotential&, const double>(),
            py::arg("collision"), py::arg("positions"),
            py::arg("barrier_potential"), py::arg("barrier_stiffness"));
}
