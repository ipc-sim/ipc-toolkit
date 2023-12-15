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
                const VertexVertexCollision&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, const double,
                const double>(),
            py::arg("collision"), py::arg("vertices"), py::arg("edges"),
            py::arg("faces"), py::arg("dhat"), py::arg("barrier_stiffness"));
}
