#include <common.hpp>

#include <ipc/collisions/tangential/edge_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_edge_vertex_tangential_collision(py::module_& m)
{
    py::class_<
        EdgeVertexTangentialCollision, EdgeVertexCandidate,
        TangentialCollision>(m, "EdgeVertexTangentialCollision")
        .def(py::init<const EdgeVertexNormalCollision&>(), py::arg("collision"))
        .def(
            py::init<
                const EdgeVertexNormalCollision&, const VectorMax12d&,
                const NormalPotential&, const double>(),
            py::arg("collision"), py::arg("positions"),
            py::arg("normal_potential"), py::arg("normal_stiffness"));
}
