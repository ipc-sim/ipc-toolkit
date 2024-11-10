#include <common.hpp>

#include <ipc/collisions/tangential/vertex_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_vertex_vertex_tangential_collision(py::module_& m)
{
    py::class_<
        VertexVertexTangentialCollision, VertexVertexCandidate,
        TangentialCollision>(m, "VertexVertexTangentialCollision")
        .def(
            py::init<const VertexVertexNormalCollision&>(),
            py::arg("collision"))
        .def(
            py::init<
                const VertexVertexNormalCollision&, const VectorMax12d&,
                const BarrierPotential&, const double>(),
            py::arg("collision"), py::arg("positions"),
            py::arg("barrier_potential"), py::arg("barrier_stiffness"));
}
