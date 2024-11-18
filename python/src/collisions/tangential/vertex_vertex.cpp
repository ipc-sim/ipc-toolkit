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
                const NormalPotential&, const double>(),
            py::arg("collision"), py::arg("positions"),
            py::arg("normal_potential"), py::arg("normal_stiffness"));
}
