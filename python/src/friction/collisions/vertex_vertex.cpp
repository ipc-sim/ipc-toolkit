#include <common.hpp>

#include <ipc/friction/collisions/vertex_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_vertex_vertex_friction_collision(py::module_& m)
{
    py::class_<
        VertexVertexFrictionCollision, VertexVertexCandidate,
        FrictionCollision>(m, "VertexVertexFrictionCollision")
        .def(py::init<const VertexVertexCollision&>(), py::arg("collision"))
        .def(
            py::init<
                const VertexVertexCollision&, const VectorMax12d&,
                const BarrierPotential&, const double>(),
            py::arg("collision"), py::arg("positions"),
            py::arg("barrier_potential"), py::arg("barrier_stiffness"));
}
